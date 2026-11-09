#!/usr/bin/env python3
"""Turn an rrnDB download into the two-column table the 16S arm reads.

Why this exists
---------------
The 16S rRNA gene sits in an *rrn* operon and bacteria carry 1-15 copies of it.
A genome with 7 copies contributes roughly 7x the 16S reads of a genome with 1
at the same cell density, so uncorrected 16S abundance overstates multi-copy
organisms in direct proportion to their copy number. The rMLST arm targets ~53
ribosomal *protein* genes, which are single copy, so comparing the two arms'
magnitudes without this correction compares two different biases.

Without a table every taxon is divided by the fallback mean of 4.2, and a
constant divisor cancels out of every ratio - it changes the absolute
``reads_assigned`` figure and nothing else. This command produces the real table.

Getting the input
-----------------
Download the pan-taxa statistics file from https://rrndb.umms.med.umich.edu/
(Downloads -> "rrnDB-5.x_pantaxa_stats_NCBI.tsv"). It has one row per taxon at
every rank with mean/stddev/min/max copy number. Only the genus and species
rows are used here; the arm resolves species first, then genus, then falls back.

Usage
-----
    bridge build-copy-numbers \\
        -i rrnDB-5.9_pantaxa_stats_NCBI.tsv \\
        -o /data/db/16S/rrnDB_copy_numbers.tsv \\
        --reference /data/db/16S/MIMt-16S_M2c_25_04_tax_bacteria_tax_in_id.fna

``--reference`` is optional but worth passing: it reports how much of MIMt the
table actually covers. A table whose names do not match the reference's is
indistinguishable, at runtime, from no table at all. Coverage is measured with
the 16S arm's own header parser.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Only these ranks are useful: the arm looks up "Genus species" then "Genus".
WANTED_RANKS = {"genus", "species"}

# rrnDB's column headers carry spaces and have changed spelling between
# releases, so they are matched by content rather than by exact name.
NAME_HINTS = ("name",)
RANK_HINTS = ("rank",)
MEAN_HINTS = ("mean",)


def _pick_column(header: list[str], hints: tuple[str, ...],
                 must_also: tuple[str, ...] = ()) -> str | None:
    """First header field containing every hint, case-insensitively."""
    for field in header:
        low = field.strip().lower()
        if all(hint in low for hint in hints) and all(m in low for m in must_also):
            return field
    return None


def read_rrndb(path: Path) -> tuple[dict[str, float], dict[str, int], set[str]]:
    """Parse rrnDB pan-taxa stats into ``{lowercased name: mean copy number}``.

    Split on tabs directly rather than through :mod:`csv`. The file has no
    quoting, and its ``sum16slist`` column holds one entry per genome in the
    taxon - for a well-sequenced species that single field runs past the csv
    module's 131,072-character limit and raises.
    """
    with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
        header_line = handle.readline().rstrip("\n").rstrip("\r")
        header = header_line.split("\t")
        if not header_line:
            raise SystemExit(f"ERROR: {path} has no header row")

        # "16S gene copy number mean" in some releases; plain "mean" in others
        # (5.10 uses the short form), so match on content rather than spelling.
        mean_col = (_pick_column(header, MEAN_HINTS, ("copy number",))
                    or _pick_column(header, MEAN_HINTS))
        name_col = _pick_column(header, NAME_HINTS)
        rank_col = _pick_column(header, RANK_HINTS)
        missing = [label for label, col in (("name", name_col), ("rank", rank_col),
                                            ("mean copy number", mean_col))
                   if col is None]
        if missing:
            raise SystemExit(
                f"ERROR: cannot find {', '.join(missing)} in {path}\n"
                f"       columns present: {', '.join(header)}")
        index = {field: position for position, field in enumerate(header)}
        rank_at, name_at, mean_at = (index[rank_col], index[name_col],
                                     index[mean_col])
        widest = max(rank_at, name_at, mean_at)

        table: dict[str, float] = {}
        genus_names: set[str] = set()
        counts = {"genus": 0, "species": 0, "skipped_rank": 0, "unparsable": 0}
        for line in handle:
            fields = line.rstrip("\n").rstrip("\r").split("\t")
            if len(fields) <= widest:
                counts["unparsable"] += 1
                continue
            rank = fields[rank_at].strip().lower()
            if rank not in WANTED_RANKS:
                counts["skipped_rank"] += 1
                continue
            name = fields[name_at].strip()
            try:
                mean = float(fields[mean_at])
            except ValueError:
                counts["unparsable"] += 1
                continue
            if not name or mean <= 0:
                counts["unparsable"] += 1
                continue
            # Lowercased: this is exactly how load_copy_numbers keys it, and how
            # aggregate._copy_number looks it up.
            table[name.lower()] = mean
            if rank == "genus":
                # Kept as a set so the coverage report can use the real
                # rank. Guessing it from the absence of a space misfiles
                # every "Candidatus X" genus as a species.
                genus_names.add(name.lower())
            counts[rank] += 1
    return table, counts, genus_names


def reference_names(path: Path) -> tuple[set[str], set[str]]:
    """Species and genus names present in the MIMt FASTA headers.

    Uses the arm's own header parser rather than a copy of it. MIMt headers are
    positional - semicolon-delimited, genus in field 7, species in field 8, with
    underscores inside the species field and an authority to strip off the genus
    - so a second implementation would disagree with the first in exactly the
    cases that matter, and the coverage figure would be fiction.
    """
    from bridge.marker16s.reference import load_reference

    species: set[str] = set()
    genera: set[str] = set()
    for taxon in load_reference(path).values():
        if taxon.genus:
            genera.add(taxon.genus.lower())
            if taxon.species_full:
                species.add(taxon.species_full.lower())
    return species, genera


def report_coverage(table: dict[str, float], table_genera: set[str],
                    reference: Path) -> None:
    """How much of the reference this table can actually put a divisor on.

    Two directions, and confusing them raises a false alarm. rrnDB only covers
    organisms with sequenced genomes - 2,280 genera in 5.10 - while MIMt covers
    4,479, many of them environmental. So "what fraction of the reference's
    genera are in the table" is capped at about 51% however perfect the
    nomenclature, so on its own it would make a healthy table look broken.

    The figure that decides what happens at runtime is **resolvable**: a
    reference taxon gets a real divisor if its own name OR its genus is in the
    table. The figure that tests *nomenclature agreement* is the reverse - how
    much of the table the reference recognises. Both are printed.

    Both are still taxon counts. What finally matters is reads-weighted, and
    that is logged per sample as ``reads_corrected_frac``: the dominant
    organisms in clinical material are exactly the well-sequenced ones, so it
    runs higher than any figure here.
    """
    species, genera = reference_names(reference)
    if not species and not genera:
        print(f"  WARNING: parsed no taxa out of {reference}; coverage unknown")
        return

    species_hit = sum(1 for name in species if name in table)
    genus_hit = sum(1 for name in genera if name in table)
    resolvable = sum(
        1 for name in species
        if name in table or name.split(" ")[0] in table)
    # Reverse direction: of the table's own entries, how many does the
    # reference know? This is what would collapse if the two databases had
    # drifted apart on naming.
    table_genus_hit = sum(1 for name in table_genera if name in genera)

    pct = lambda hit, total: 100 * hit / max(total, 1)  # noqa: E731

    print("\nCoverage against the reference")
    print(f"  resolvable at species or genus : {resolvable:>6} / {len(species):<6} "
          f"({pct(resolvable, len(species)):.1f}%)   <- what the arm achieves")
    print(f"  exact species matches          : {species_hit:>6} / {len(species):<6} "
          f"({pct(species_hit, len(species)):.1f}%)")
    print(f"  reference genera in the table  : {genus_hit:>6} / {len(genera):<6} "
          f"({pct(genus_hit, len(genera)):.1f}%)   (capped by rrnDB's size)")
    print(f"  table genera in the reference  : {table_genus_hit:>6} / "
          f"{len(table_genera):<6} ({pct(table_genus_hit, len(table_genera)):.1f}%)"
          f"   <- nomenclature agreement")

    if len(species) and resolvable / len(species) < 0.5:
        print("\n  WARNING: under half the reference's species can be given a real")
        print("           copy number. Most taxa will fall back to the 4.2 default,")
        print("           which cancels out of every ratio and corrects nothing.")
    elif len(table_genera) and table_genus_hit / len(table_genera) < 0.6:
        print("\n  WARNING: the reference recognises under 60% of the table's own")
        print("           genera, which suggests the rrnDB release and the")
        print("           reference disagree on nomenclature rather than simply")
        print("           covering different organisms.")
    else:
        print("\n  Coverage is sound. The per-sample logs will report")
        print("  reads_corrected_frac, which is the figure that matters.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="bridge build-copy-numbers",
        description="Build the 16S rrn copy-number table from an rrnDB download.")
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="rrnDB-5.x_pantaxa_stats_NCBI.tsv")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="where to write name<TAB>copy_number "
                             "(point references.s16_copy_numbers at it)")
    parser.add_argument("--reference", type=Path, default=None,
                        help="MIMt FASTA; reports how much of it the table covers")
    args = parser.parse_args(argv)

    if not args.input.is_file():
        print(f"ERROR: no such file: {args.input}", file=sys.stderr)
        return 1

    table, counts, table_genera = read_rrndb(args.input)
    if not table:
        print(f"ERROR: no genus or species rows found in {args.input}",
              file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        handle.write("# rrn 16S copy numbers, built from "
                     f"{args.input.name} by bridge build-copy-numbers\n")
        handle.write("# name\tcopy_number\n")
        for name in sorted(table):
            handle.write(f"{name}\t{table[name]:.4g}\n")

    values = sorted(table.values())
    print(f"Wrote {args.output}")
    print(f"  species rows : {counts['species']:>6}")
    print(f"  genus rows   : {counts['genus']:>6}")
    print(f"  other ranks skipped: {counts['skipped_rank']:>6}")
    if counts["unparsable"]:
        print(f"  unparsable rows    : {counts['unparsable']:>6}")
    print(f"  copy number range  : {values[0]:.2f} - {values[-1]:.2f} "
          f"(median {values[len(values) // 2]:.2f})")

    if args.reference:
        if args.reference.is_file():
            report_coverage(table, table_genera, args.reference)
        else:
            print(f"  WARNING: no such reference: {args.reference}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
