"""The ResFinder database: template names, gene families, phenotypes, the KMA index.

    python -m bridge.amr.cli build-db --db-dir /data/db/resfinder_db

Template names
--------------
ResFinder names every sequence ``<gene>_<variant>_<accession>``, e.g.
``blaCTX-M-15_1_AY044436``. The gene is what a report names; the variant and
accession tell apart submissions of the same allele. The gene is taken up to
the first ``_<digits>_`` because gene names themselves never contain one
while some accessions do (``NG_047471``).

Families
--------
A *family* is the level at which two questions are asked that an allele is
too fine for: is this gene intrinsic to the organism it was found with
(``host.py``), and does the allele matter for the phenotype (the spectrum
check below). ``family_of`` strips the allele designator by rule, after an
explicit override table (``data/family_overrides.tsv``) that holds the
exceptions - the same enzyme under two naming conventions, mosaic tet genes,
and families that must NOT be merged (fosA7 is Salmonella's chromosomal gene;
fosA is Klebsiella's and Enterobacter's).

The rule is deliberately conservative. A name it does not recognise is its
own family, which costs nothing but an ``undetermined`` host verdict; a rule
that merged too eagerly would hand one family's verdict to another.

Phenotypes
----------
``phenotypes.txt`` gives the antibiotics each allele confers resistance to,
and it is allele-specific where it matters: blaSHV-1 is a narrow-spectrum
penicillinase, blaSHV-12 an ESBL. So a call that did not resolve the allele
can only borrow its nearest allele's list if the family agrees with itself.
:meth:`PhenotypeDB.family_spectrum` answers that - the union and intersection
of the antibiotic lists across every allele of the family.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import re
import subprocess
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

LOGGER = logging.getLogger("bridge.amr.database")

TEMPLATE_RE = re.compile(r"^(?P<gene>.+?)_(?P<variant>\d+)_(?P<accession>.+)$")


@dataclass(frozen=True)
class TemplateName:
    template: str
    gene: str
    variant: str
    accession: str


def parse_template(name: str) -> TemplateName:
    name = str(name).strip()
    match = TEMPLATE_RE.match(name)
    if not match:
        return TemplateName(name, name, "", "")
    return TemplateName(name, match["gene"], match["variant"], match["accession"])


# ---------------------------------------------------------------------------
# Families
# ---------------------------------------------------------------------------

# Ordered: the first matching rule wins. Each maps a gene name to its family.
_FAMILY_RULES: tuple[tuple[re.Pattern, str], ...] = tuple(
    (re.compile(pattern), repl) for pattern, repl in (
        # beta-lactamases: strip the allele number ('blaCTX-M-15' -> 'blaCTX-M',
        # 'blaTEM-1B' -> 'blaTEM', 'blaOXA-48' -> 'blaOXA', 'blaOXY-1-1' -> 'blaOXY')
        (r"^(bla[A-Za-z0-9]+?(?:-[A-Z]+)?)-\d+[A-Za-z]?(?:-\d+)?$", r"\1"),
        # modifying enzymes: the enzyme class and site ('aac(6')-Ib-cr' -> "aac(6')")
        (r"^((?:aac|aph|ant|apm)\([^)]*\)).*$", r"\1"),
        (r"^(aad[A-Z])\d*$", r"\1"),
        (r"^(qnr[A-Z])\d+$", r"\1"),
        (r"^(dfr[A-Z])\d*[a-z]?$", r"\1"),
        (r"^(cat[AB])\d+$", r"\1"),
        (r"^(cml[AB])\d*$", r"\1"),
        (r"^(mcr)-[\d.]+$", r"\1"),
        (r"^(tet\(X)\d*\)$", r"\1)"),
        (r"^(tmex[CD]|toprJ|qepA|oqx[AB])\d*$", r"\1"),
        (r"^(ARR)-\d+$", r"\1"),
    )
)


def load_family_overrides(path: Path | None) -> list[tuple[re.Pattern, str]]:
    """``gene_regex<TAB>family`` rows; '#' comments and blank lines skipped."""
    rows: list[tuple[re.Pattern, str]] = []
    if path is None or not Path(path).is_file():
        return rows
    header_seen = False
    with open(path, encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cells = [cell.strip() for cell in line.split("\t")]
            if not header_seen:
                header_seen = True
                if cells[:2] != ["gene_regex", "family"]:
                    raise ValueError(f"{Path(path).name}:{lineno}: header must start "
                                     f"gene_regex<TAB>family, got {cells[:2]}")
                continue
            if len(cells) < 2 or not cells[0] or not cells[1]:
                raise ValueError(f"{Path(path).name}:{lineno}: need gene_regex and family")
            rows.append((re.compile(cells[0]), cells[1]))
    return rows


def family_of(gene: str, overrides: list[tuple[re.Pattern, str]] | None = None) -> str:
    """Gene name -> family key. Overrides first (full match), then the rules."""
    gene = str(gene).strip()
    for pattern, family in overrides or ():
        if pattern.fullmatch(gene):
            return family
    for pattern, repl in _FAMILY_RULES:
        if pattern.match(gene):
            return pattern.sub(repl, gene)
    return gene


# ---------------------------------------------------------------------------
# Phenotypes
# ---------------------------------------------------------------------------

def normalise_agent(name: str) -> str:
    """Antibiotic name -> comparison key.

    ResFinder writes 'Amoxicillin+Clavulanic acid', EUCAST
    'amoxicillin-clavulanic acid'; both become 'amoxicillin-clavulanic acid'.
    """
    text = str(name).strip().lower()
    text = re.sub(r"\s*(\+|/|\band\b)\s*", "-", text)
    text = re.sub(r"\s*-\s*", "-", text)
    return re.sub(r"\s+", " ", text)


@dataclass(frozen=True)
class PhenotypeEntry:
    template: str
    gene: str
    drug_class: str
    antibiotics: tuple[str, ...]
    mechanism: str = ""
    notes: str = ""
    pmid: str = ""
    required_gene: str = ""


def _column(header: list[str], *needles: str) -> int | None:
    lowered = [cell.strip().lower() for cell in header]
    for needle in needles:
        for index, cell in enumerate(lowered):
            if needle in cell:
                return index
    return None


def _split_agents(text: str) -> tuple[str, ...]:
    agents = []
    for part in str(text).split(","):
        part = part.strip()
        if part and part.lower() not in ("unknown", "none", "-"):
            agents.append(part)
    return tuple(agents)


@dataclass
class PhenotypeDB:
    by_template: dict[str, PhenotypeEntry]
    overrides: list[tuple[re.Pattern, str]] = field(default_factory=list)
    _by_gene: dict[str, list[PhenotypeEntry]] = field(default_factory=dict, init=False)
    _spectra: dict[str, tuple[frozenset, frozenset]] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        for entry in self.by_template.values():
            self._by_gene.setdefault(entry.gene, []).append(entry)

    def __len__(self) -> int:
        return len(self.by_template)

    def lookup(self, template: str) -> PhenotypeEntry | None:
        """This template's entry; else any entry for the same gene."""
        template = str(template).strip()
        if template in self.by_template:
            return self.by_template[template]
        entries = self._by_gene.get(parse_template(template).gene)
        return entries[0] if entries else None

    def family_spectrum(self, family: str) -> tuple[frozenset, frozenset]:
        """(union, intersection) of normalised antibiotics across the family."""
        if family not in self._spectra:
            sets = [frozenset(normalise_agent(a) for a in entry.antibiotics)
                    for gene, entries in self._by_gene.items()
                    if family_of(gene, self.overrides) == family
                    for entry in entries]
            if sets:
                union = frozenset().union(*sets)
                inter = frozenset.intersection(*sets)
            else:
                union = inter = frozenset()
            self._spectra[family] = (union, inter)
        return self._spectra[family]


def load_phenotypes(path: Path, overrides: list[tuple[re.Pattern, str]] | None = None) -> PhenotypeDB:
    """ResFinder ``phenotypes.txt`` (tab-separated, one header row)."""
    entries: dict[str, PhenotypeEntry] = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = None
        for row in reader:
            if not row or not any(cell.strip() for cell in row):
                continue
            if header is None:
                header = row
                cols = {
                    "template": _column(header, "gene_accession", "gene accession") or 0,
                    "class": _column(header, "class"),
                    "phenotype": _column(header, "phenotype"),
                    "pmid": _column(header, "pmid"),
                    "mechanism": _column(header, "mechanism"),
                    "notes": _column(header, "notes"),
                    "required": _column(header, "required"),
                }
                if cols["phenotype"] is None:
                    raise ValueError(f"{path}: no Phenotype column in header {header}")
                continue

            def cell(key: str) -> str:
                index = cols[key]
                return row[index].strip() if index is not None and index < len(row) else ""

            template = cell("template")
            if not template:
                continue
            entries.setdefault(template, PhenotypeEntry(
                template=template,
                gene=parse_template(template).gene,
                drug_class=cell("class"),
                antibiotics=_split_agents(cell("phenotype")),
                mechanism=cell("mechanism"),
                notes=cell("notes"),
                pmid=cell("pmid"),
                required_gene=cell("required"),
            ))
    if not entries:
        raise ValueError(f"{path}: no phenotype rows read")
    return PhenotypeDB(entries, list(overrides or []))


# ---------------------------------------------------------------------------
# Building the KMA index
# ---------------------------------------------------------------------------

# Files in a resfinder_db checkout that are not class FASTAs.
_SKIP_FASTA = {"all.fsa", "resfinder.fsa"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _git(db_dir: Path, *args: str) -> str | None:
    try:
        result = subprocess.run(["git", "-C", str(db_dir), *args], capture_output=True,
                                text=True, check=False)
    except OSError:
        return None
    return result.stdout.strip() or None if result.returncode == 0 else None


def concatenate_fastas(db_dir: Path, out_fasta: Path) -> tuple[list[Path], int]:
    """Every class FASTA into one, first copy of each template name kept.

    ResFinder ships one FASTA per antibiotic class. A template in two of them
    is the same sequence filed twice; KMA would index it twice and split its
    reads between the copies, so only the first is kept (and a copy that
    differs in sequence is logged, because that would be a database fault).
    """
    sources = sorted(path for path in Path(db_dir).glob("*.fsa")
                     if path.name not in _SKIP_FASTA)
    if not sources:
        raise FileNotFoundError(f"no *.fsa class files in {db_dir}")
    seen: dict[str, str] = {}
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(out_fasta, "w", encoding="utf-8", newline="\n") as out:
        for source in sources:
            name, chunks = None, []

            def flush() -> None:
                nonlocal count
                if name is None:
                    return
                sequence = "".join(chunks).upper()
                if name in seen:
                    if seen[name] != sequence:
                        LOGGER.warning("template %s appears twice with different "
                                       "sequences; the first (%s) is kept", name, source.name)
                    return
                seen[name] = sequence
                out.write(f">{name}\n{sequence}\n")
                count += 1

            with open(source, encoding="utf-8", errors="replace") as handle:
                for raw in handle:
                    line = raw.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        flush()
                        name, chunks = line[1:].split()[0], []
                    else:
                        chunks.append(line)
                flush()
    return sources, count


def kma_index_command(kma: str, kma_index: str | None = None) -> list[str]:
    """``kma_index`` beside ``kma`` if it is there, else ``kma index``."""
    from bridge.tools import companion_executable
    if kma_index:
        return [kma_index]
    try:
        return [companion_executable("kma_index", kma)]
    except (FileNotFoundError, OSError):
        return [kma, "index"]


def build_db(db_dir: Path, prefix: Path | None = None, *, kma: str = "kma",
             kma_index: str | None = None, dry_run: bool = False) -> dict:
    """Concatenate, index and record a ResFinder checkout. Returns the manifest."""
    from .settings import DB_MANIFEST_NAME

    db_dir = Path(db_dir).resolve()
    prefix = Path(prefix) if prefix else db_dir / "kma" / "resfinder"
    phenotypes = db_dir / "phenotypes.txt"
    if not phenotypes.is_file():
        raise FileNotFoundError(f"{phenotypes} not found: is {db_dir} a resfinder_db checkout?")
    fasta = prefix.with_suffix(".fsa")
    sources, n_templates = concatenate_fastas(db_dir, fasta)
    command = kma_index_command(kma, kma_index) + ["-i", str(fasta), "-o", str(prefix)]
    if not dry_run:
        from bridge.tools import run
        run(command, "kma_index")

    manifest = {
        "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "db_dir": str(db_dir),
        "index_prefix": str(prefix),
        "git_commit": _git(db_dir, "rev-parse", "HEAD"),
        "git_date": _git(db_dir, "log", "-1", "--format=%cI"),
        "version_file": ((db_dir / "VERSION").read_text().strip()
                         if (db_dir / "VERSION").is_file() else None),
        "templates": n_templates,
        "fasta_sha256": _sha256(fasta),
        "phenotypes": str(phenotypes),
        "phenotypes_sha256": _sha256(phenotypes),
        "sources": {path.name: _sha256(path) for path in sources},
        "command": command,
    }
    (prefix.parent / DB_MANIFEST_NAME).write_text(json.dumps(manifest, indent=2) + "\n",
                                                  encoding="utf-8")
    return manifest


def read_manifest(path: Path | None) -> dict | None:
    if path is None or not Path(path).is_file():
        return None
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
