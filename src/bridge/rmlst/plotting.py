from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


EVIDENCE_COLUMNS = {
    "allele": (
        ("allele_evidence_unambiguous_loci", "unambiguous_loci"),
        ("allele_evidence_ambiguous_loci", "ambiguous_loci"),
    ),
    "em": (
        ("em_evidence_unambiguous_loci", "unambiguous_loci"),
        ("em_evidence_ambiguous_loci", "ambiguous_loci"),
    ),
}


def _first_existing(columns: pd.Index, candidates: tuple[str, ...]) -> str:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise ValueError(f"None of the required columns were found: {', '.join(candidates)}")


def _reported_mask(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _short_species_name(name: str) -> str:
    parts = str(name).split()
    if len(parts) >= 2:
        return f"{parts[0][0]}. {' '.join(parts[1:])}"
    return str(name)


def prepare_called_species(
    calls: pd.DataFrame,
    evidence_source: str = "allele",
    log_floor: float = 1e-6,
) -> pd.DataFrame:
    if evidence_source not in EVIDENCE_COLUMNS:
        raise ValueError(f"Unknown evidence source: {evidence_source}")
    if log_floor <= 0:
        raise ValueError("log_floor must be greater than zero")
    required = {"species", "reported", "pi"}
    missing = sorted(required - set(calls.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    called = calls.loc[_reported_mask(calls["reported"])].copy()
    if called.empty:
        return called

    unambiguous_candidates, ambiguous_candidates = EVIDENCE_COLUMNS[evidence_source]
    unambiguous_column = _first_existing(calls.columns, unambiguous_candidates)
    ambiguous_column = _first_existing(calls.columns, ambiguous_candidates)
    called["pi"] = pd.to_numeric(called["pi"], errors="raise")
    if (called["pi"] < 0).any():
        raise ValueError("pi contains negative values, which cannot use a log axis")
    # A species called from its allele profile that never reached EM has no
    # abundance estimate, and carries NaN rather than 0.0 to say so. Draw it at
    # the floor like a zero: leaving it NaN makes matplotlib skip the bar
    # entirely, silently dropping a reported organism from the figure.
    called["plot_pi"] = called["pi"].fillna(0.0).clip(lower=log_floor)
    called["unambiguous_loci"] = (
        pd.to_numeric(called[unambiguous_column], errors="coerce").fillna(0).astype(int)
    )
    called["ambiguous_loci"] = (
        pd.to_numeric(called[ambiguous_column], errors="coerce").fillna(0).astype(int)
    )
    called["species_label"] = called["species"].map(_short_species_name)
    return called


def default_plot_title(sample_id: str) -> str:
    readable = str(sample_id).replace("__", " | ").replace("_", " ")
    return "Called species abundance\n" + textwrap.fill(readable, width=68)


def write_called_species_pi_pdf(
    calls: pd.DataFrame,
    output_path: Path,
    sample_id: str,
    *,
    title: str | None = None,
    evidence_source: str = "allele",
    sort_order: str = "pi",
    log_floor: float = 1e-6,
) -> Path:
    if sort_order not in {"pi", "species"}:
        raise ValueError("sort_order must be pi or species")
    called = prepare_called_species(calls, evidence_source, log_floor)
    if sort_order == "pi" and not called.empty:
        called = called.sort_values(["pi", "species"], ascending=[False, True])
    elif not called.empty:
        called = called.sort_values("species")
    called = called.reset_index(drop=True)
    plot_title = textwrap.fill(title, width=68) if title else default_plot_title(sample_id)

    sns.set_theme(style="whitegrid", context="notebook", font_scale=1.0)
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.size": 12,
        "axes.titlesize": 18,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "text.usetex": False,
    })

    figure_width = max(8.0, min(16.0, 2.0 + 1.65 * max(len(called), 1)))
    fig, ax = plt.subplots(figsize=(figure_width, 6.0))
    fig.subplots_adjust(left=0.14, right=0.98, top=0.80, bottom=0.27)
    ax.set_title(plot_title, loc="left", fontweight="bold", pad=14)

    if called.empty:
        ax.set_axis_off()
        ax.text(
            0.5, 0.48, "No species were reported.",
            transform=ax.transAxes, ha="center", va="center", fontsize=16,
            color="#5E6872",
        )
    else:
        order = called["species_label"].tolist()
        palette = sns.color_palette("colorblind", n_colors=len(called))
        sns.barplot(
            data=called, x="species_label", y="plot_pi", order=order,
            hue="species_label", hue_order=order, palette=palette, dodge=False,
            errorbar=None, edgecolor="#30363d", linewidth=0.9, legend=False, ax=ax,
        )
        positive = called.loc[called["pi"] > 0, "pi"]
        minimum = max(
            log_floor * 0.5,
            positive.min() / 3 if not positive.empty else log_floor / 2,
        )
        maximum = max(1.2, called["plot_pi"].max() * 3.0)
        ax.set_yscale("log")
        ax.set_ylim(minimum, maximum)
        ax.set_xlabel("Called species")
        ax.set_ylabel("Estimated abundance, pi (log scale)")
        ax.grid(axis="x", visible=False)
        ax.grid(axis="y", which="major", color="#D9DEE3", linewidth=0.8)
        ax.grid(axis="y", which="minor", color="#EEF1F3", linewidth=0.45)
        sns.despine(ax=ax)
        if len(called) > 4 or any(len(label) > 20 for label in order):
            plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")
        for patch, row in zip(ax.patches, called.itertuples(index=False)):
            x = patch.get_x() + patch.get_width() / 2
            ax.text(
                x, row.plot_pi * 1.22,
                f"U: {row.unambiguous_loci}\nA: {row.ambiguous_loci}",
                ha="center", va="bottom", fontsize=12, color="#20252B",
            )
        ax.text(
            1.0, -0.22,
            f"Labels: U = unambiguous loci; A = ambiguous loci ({evidence_source} evidence).",
            transform=ax.transAxes, ha="right", va="top", fontsize=12,
            color="#5E6872",
        )
        floored = []
        if (called["pi"] == 0).any():
            floored.append("Zero")
        if called["pi"].isna().any():
            # Distinct from zero, and the distinction is the point: this
            # organism's abundance was never estimated, not measured as absent.
            floored.append("Unestimated")
        if floored:
            ax.text(
                0.0, -0.22,
                f"{' and '.join(floored).lower().capitalize()} pi values "
                f"displayed at {log_floor:g}.",
                transform=ax.transAxes, ha="left", va="top", fontsize=12,
                color="#5E6872",
            )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="pdf", metadata={"Title": plot_title})
    plt.close(fig)
    return output_path
