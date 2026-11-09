"""Shared synthetic inputs for the AMR arm tests: a tiny ResFinder, KMA output, a results tree."""

from __future__ import annotations

import csv
import json
from pathlib import Path

PHENOTYPES = """\
Gene_accession no.\tClass\tPhenotype\tPMID\tMechanism of resistance\tNotes\tRequired_gene
blaTEM-1B_1_AY458016\tBeta-lactam\tAmoxicillin, Ampicillin, Cephalothin, Piperacillin, Ticarcillin\t123\tAntibiotic inactivation\t\t
blaTEM-1A_1_HM749966\tBeta-lactam\tAmoxicillin, Ampicillin, Cephalothin, Piperacillin, Ticarcillin\t123\tAntibiotic inactivation\t\t
blaSHV-1_1_AF148850\tBeta-lactam\tAmoxicillin, Ampicillin, Piperacillin, Ticarcillin\t1\tAntibiotic inactivation\t\t
blaSHV-12_1_AF462395\tBeta-lactam\tAmoxicillin, Ampicillin, Aztreonam, Cefotaxime, Ceftazidime, Ceftriaxone, Piperacillin, Ticarcillin\t2\tAntibiotic inactivation\tESBL\t
blaCTX-M-15_1_AY044436\tBeta-lactam\tAmoxicillin, Ampicillin, Aztreonam, Cefepime, Cefotaxime, Ceftazidime, Ceftriaxone, Piperacillin, Ticarcillin\t3\tAntibiotic inactivation\t\t
blaACT-1_1_U58495\tBeta-lactam\tAmoxicillin, Amoxicillin+Clavulanic acid, Ampicillin, Cefoxitin\t4\tAntibiotic inactivation\t\t
tet(M)_1_X92947\tTetracycline\tDoxycycline, Minocycline, Tetracycline\t5\tTarget protection\t\t
"""

RES_HEADER = ("#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\t"
              "Query_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n")
MAPSTAT_HEADER = ("## method\tKMA\n## version\t1.6.6\n## database\tresfinder\n"
                  "## fragmentCount\t1000000\n"
                  "# refSequence\treadCount\tfragmentCount\tmapScoreSum\trefCoveredPositions\t"
                  "refConsensusSum\tbpTotal\tdepthVariance\tnucHighDepthVariance\tdepthMax\t"
                  "snpSum\tinsertSum\tdeletionSum\treadCountAln\tfragmentCountAln\n")


def res_line(template, identity=100.0, coverage=100.0, depth=50.0, length=861,
             query_identity=None):
    qi = identity if query_identity is None else query_identity
    # KMA pads the template column with spaces; the parser must cope.
    return (f"{template:<40}\t   44521\t      12\t{length:>8}\t{identity:7.2f}\t{coverage:7.2f}\t"
            f"{qi:7.2f}\t{coverage:7.2f}\t{depth:8.2f}\t44456.38\t1.0e-26\n")


def mapstat_line(template, fragments, reads=None):
    reads = fragments * 2 if reads is None else reads
    return f"{template}\t{reads}\t{fragments}\t9999\t861\t861\t43000\t1.0\t0\t80\t0\t0\t0\t{reads}\t{fragments}\n"


def write_phenotypes(folder: Path) -> Path:
    path = Path(folder) / "phenotypes.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(PHENOTYPES, encoding="utf-8")
    return path


def write_kma_output(prefix: Path, hits: list[tuple]) -> None:
    """hits: (template, identity, coverage, depth, fragments[, query_identity])."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    res = RES_HEADER
    mapstat = MAPSTAT_HEADER
    for hit in hits:
        template, identity, coverage, depth, fragments = hit[:5]
        qi = hit[5] if len(hit) > 5 else None
        res += res_line(template, identity, coverage, depth, query_identity=qi)
        mapstat += mapstat_line(template, fragments)
    Path(f"{prefix}.res").write_text(res, encoding="utf-8")
    Path(f"{prefix}.mapstat").write_text(mapstat, encoding="utf-8")


def amr_sample(root: Path, name: str, hits: list[tuple], total: int, settings) -> Path:
    """An arm result as the arm writes it, without running KMA."""
    from bridge.amr.calls import classify
    from bridge.amr.database import load_family_overrides, load_phenotypes
    from bridge.amr.kma import merge_hits, parse_mapstat, parse_res, write_hits
    from bridge.amr.phenotype import annotate

    amr_dir = Path(root) / name / "amr"
    write_kma_output(amr_dir / "kma", hits)
    overrides = load_family_overrides(settings.family_overrides)
    db = load_phenotypes(settings.phenotypes, overrides)
    rows = merge_hits(name, parse_res(amr_dir / "kma.res"),
                      parse_mapstat(amr_dir / "kma.mapstat")[1], overrides)
    classify(rows, settings.thresholds)
    for row in rows:
        annotate(row, db)
    write_hits(rows, amr_dir / "amr_hits.csv")
    (amr_dir / "amr_manifest.json").write_text(json.dumps({"total_read_pairs": total}))
    return amr_dir


def id_summary(root: Path, name: str, organisms: list[tuple[str, str, float]]) -> None:
    """(organism, basis, median_projected_depth) rows, as reconcile writes them."""
    folder = Path(root) / name / "bacterial_ID"
    folder.mkdir(parents=True, exist_ok=True)
    columns = ["sample", "run_id", "organism", "genus", "basis", "status",
               "median_projected_depth"]
    with open(folder / "identification_summary.csv", "w", newline="", encoding="utf-8") as h:
        writer = csv.DictWriter(h, fieldnames=columns)
        writer.writeheader()
        for organism, basis, depth in organisms:
            writer.writerow({"sample": name, "run_id": "RUN1", "organism": organism,
                             "genus": organism.split()[0], "basis": basis,
                             "status": "report", "median_projected_depth": depth})


def settings_with_db(tmp: Path):
    """The repository's AMR settings, pointed at the synthetic phenotypes.txt."""
    from dataclasses import replace
    from bridge.amr.settings import load_amr_settings
    phenotypes = write_phenotypes(Path(tmp) / "resfinder_db")
    return replace(load_amr_settings(), phenotypes=phenotypes,
                   db=Path(tmp) / "resfinder_db" / "kma" / "resfinder")


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
