"""Create summary plots from an angiogenesis monitoring CSV.

Example:
    python assignment3/scripts/plot_angiogenesis_metrics.py \
        --input assignment3/Angiogenesis/results/angiogenesis_metrics.csv

The script reads a CSV produced by `MonitoringSteppable` and saves grouped plots
into a `plots/` directory next to the CSV by default.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd

SMOOTH_WINDOW = 15  # rolling-average window for noisy / rate series
ZERO_CROSSING_MIN_GAP_FRACTION = 0.07


DEFAULT_INPUT = Path("../Angiogenesis/results/angiogenesis_metrics.csv")
DEFAULT_FORMAT = "png"
DEFAULT_DPI = 160


def _smooth(series: pd.Series, window: int = SMOOTH_WINDOW) -> pd.Series:
    """Rolling mean with min_periods=1 so edges do not go NaN."""
    return series.rolling(window, center=True, min_periods=1).mean()


def _twin_right(ax: plt.Axes, color: str) -> plt.Axes:
    """Return a right-hand twin axis styled with *color*."""
    ax2 = ax.twinx()
    ax2.spines["right"].set_color(color)
    ax2.tick_params(axis="y", colors=color)
    ax2.yaxis.label.set_color(color)
    return ax2


def _decorate(
    ax: plt.Axes,
    *,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
) -> None:
    if title:
        ax.set_title(title, fontsize=9)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25, linestyle="--")


def _debounce_crossings(crossings: list[int], min_gap: int) -> list[int]:
    if not crossings:
        return []
    min_gap = max(1, int(min_gap))
    kept = [int(crossings[0])]
    for cx in crossings[1:]:
        if int(cx) - kept[-1] >= min_gap:
            kept.append(int(cx))
    return kept


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate plots from an angiogenesis monitoring CSV.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Path to the angiogenesis metrics CSV file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for generated figures. Defaults to <csv_dir>/plots.",
    )
    parser.add_argument(
        "--prefix",
        default="angiogenesis",
        help="Filename prefix for generated figures.",
    )
    parser.add_argument(
        "--format",
        default=DEFAULT_FORMAT,
        choices=("png", "pdf", "svg"),
        help="Figure format.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help="Figure DPI for raster output.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display figures interactively in addition to saving them.",
    )
    return parser.parse_args()


def load_metrics(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df = pd.read_csv(csv_path)
    if df.empty:
        raise ValueError(f"CSV file is empty: {csv_path}")
    if "mcs" not in df.columns:
        raise ValueError("CSV must contain an 'mcs' column.")

    df = df.sort_values("mcs").reset_index(drop=True)
    return df


def resolve_output_dir(csv_path: Path, output_dir: Path | None) -> Path:
    resolved = output_dir if output_dir is not None else csv_path.parent / "plots"
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def available_columns(df: pd.DataFrame, columns: Iterable[str]) -> list[str]:
    return [column for column in columns if column in df.columns]


def save_figure(fig: plt.Figure, output_dir: Path, prefix: str, stem: str, fmt: str, dpi: int) -> Path:
    output_path = output_dir / f"{prefix}_{stem}.{fmt}"
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output_path


# ── 1. Tumor dynamics ─────────────────────────────────────────────────────────

def plot_tumor_dynamics(df: pd.DataFrame) -> plt.Figure | None:
    """2x2 grid: cell counts · phenotype fractions · volumes · growth pressure."""
    has_counts = any(c in df.columns for c in ["tumor_like_cells", "normal_cells", "hypoxic_cells", "necrotic_cells"])
    has_fractions = any(c in df.columns for c in ["hypoxic_fraction", "necrotic_fraction"])
    has_volumes = any(c in df.columns for c in ["total_tumor_volume", "avg_tumor_volume"])
    has_pressure = "avg_tumor_target_volume" in df.columns and "avg_tumor_volume" in df.columns
    if not any([has_counts, has_fractions, has_volumes, has_pressure]):
        return None

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    fig.suptitle("1 · Tumor Dynamics", fontsize=12, fontweight="bold")
    mcs = df["mcs"]

    # ── counts ────────────────────────────────────────────────────────────────
    ax = axes[0, 0]
    CELL_STYLES: list[tuple[str, str, str, float]] = [
        ("tumor_like_cells", "Total tumor",  "black",       2.0),
        ("normal_cells",     "Normal",        "tab:blue",    1.2),
        ("hypoxic_cells",    "Hypoxic",       "tab:orange",  1.2),
        ("necrotic_cells",   "Necrotic",      "tab:red",     1.2),
    ]
    plotted = False
    for col, label, color, lw in CELL_STYLES:
        if col in df.columns:
            ax.plot(mcs, df[col], label=label, color=color, linewidth=lw)
            plotted = True
    _decorate(ax,
              title="Cell counts\n(slow growth → acceleration → O\u2082 limitation)",
              ylabel="Number of cells")
    if plotted:
        ax.legend(fontsize=8)

    # ── phenotype fractions + VEGF2 right axis ────────────────────────────────
    ax = axes[0, 1]
    frac_plotted = False
    for col, label, color in [
        ("hypoxic_fraction",  "Hypoxic fraction",  "tab:orange"),
        ("necrotic_fraction", "Necrotic fraction", "tab:red"),
    ]:
        if col in df.columns:
            ax.plot(mcs, df[col], label=label, color=color, linewidth=1.4)
            frac_plotted = True
    _decorate(ax,
              title="Phenotype fractions & VEGF2\n(hypoxia drives VEGF2 → angiogenesis → O\u2082 recovery → fractions drop)",
              ylabel="Fraction")
    if frac_plotted:
        ax.set_ylim(bottom=0)
        ax.legend(loc="upper left", fontsize=8)
    if "mean_tumor_vegf2" in df.columns:
        ax2 = _twin_right(ax, "tab:purple")
        ax2.plot(mcs, _smooth(df["mean_tumor_vegf2"]),
                 color="tab:purple", linestyle="--", linewidth=1.2, label="Mean VEGF2 (smooth)")
        ax2.set_ylabel("Mean VEGF2")
        ax2.legend(loc="upper right", fontsize=8)

    # ── volumes ───────────────────────────────────────────────────────────────
    ax = axes[1, 0]
    if "total_tumor_volume" in df.columns:
        ax.plot(mcs, df["total_tumor_volume"], color="tab:blue",
                label="Total volume", linewidth=1.5)
        ax.set_ylabel("Total volume (voxels)")
    if "avg_tumor_volume" in df.columns:
        ax2 = _twin_right(ax, "tab:cyan")
        ax2.plot(mcs, df["avg_tumor_volume"], color="tab:cyan",
                linestyle="--", linewidth=1.2, label="Avg cell volume")
        ax2.set_ylabel("Avg cell volume (voxels)")
        ax2.legend(loc="upper left", fontsize=8)
    
    # ── growth pressure: target vs actual ─────────────────────────────────────
    ax = axes[1, 1]
    pres_plotted = False
    if has_pressure:
        ax.plot(mcs, df["avg_tumor_target_volume"], color="tab:green",
                linewidth=1.5, label="Avg target volume")
        ax.plot(mcs, df["avg_tumor_volume"], color="tab:blue",
                linestyle="--", linewidth=1.2, label="Avg actual volume")
        ax.fill_between(mcs, df["avg_tumor_volume"], df["avg_tumor_target_volume"],
                        alpha=0.15, color="tab:green", label="Growth pressure")
        pres_plotted = True
    _decorate(ax, title="Growth pressure (target \u2212 actual)\n+ HIF-1\u03b1 if network enabled",
              xlabel="MCS", ylabel="Volume (voxels)")
    if pres_plotted:
        ax.legend(loc="upper left", fontsize=8)
    # overlay HIF-1a if present
    if "mean_tumor_hif1a" in df.columns:
        ax2 = _twin_right(ax, "tab:purple")
        ax2.plot(mcs, df["mean_tumor_hif1a"], color="tab:purple",
                 linewidth=1.2, linestyle=":", label="Mean HIF-1\u03b1")
        ax2.set_ylabel("Mean HIF-1\u03b1")
        ax2.legend(loc="lower right", fontsize=8)

    fig.tight_layout()
    return fig


# ── 2. Vascular response ──────────────────────────────────────────────────────

def plot_vascular_response(df: pd.DataFrame) -> plt.Figure | None:
    """Panels: endothelial populations · cell sizes · endothelial growth rate."""
    # Only include endothelial-related columns
    pop_cols = available_columns(
        df, ["active_neovascular_cells", "inactive_neovascular_cells"]
    )
    size_cols = available_columns(df, ["avg_endothelial_volume"])
    rate_cols = available_columns(df, ["avg_endothelial_volume_growth_rate"])
    if not pop_cols and not size_cols and not rate_cols:
        return None

    n_panels = sum([bool(pop_cols), bool(size_cols), bool(rate_cols)])
    fig, axes = plt.subplots(n_panels, 1, figsize=(11, 3.5 * n_panels), sharex=True)
    axes = [axes] if n_panels == 1 else list(axes)
    fig.suptitle("2 · Endothelial Response", fontsize=12, fontweight="bold")
    mcs = df["mcs"]
    panel = 0

    # populations
    if pop_cols:
        ax = axes[panel]; panel += 1
        ENDO_STYLES: list[tuple[str, str, str, float]] = [
            ("active_neovascular_cells",   "Active (tip cell)", "tab:green",  1.3),
            ("inactive_neovascular_cells", "Inactive (stalk)",  "tab:olive",  1.3),
        ]
        for col, label, color, lw in ENDO_STYLES:
            if col in df.columns:
                ax.plot(mcs, df[col], label=label, color=color, linewidth=lw)
        _decorate(ax,
                  title="Endothelial populations\n(VEGF2 chemotaxis activates sprouting; active count rises with hypoxia)",
                  ylabel="Cell count")
        ax.legend(fontsize=8)

    # mean cell volume
    if size_cols:
        ax = axes[panel]; panel += 1
        for col in size_cols:
            ax.plot(mcs, df[col], label=col.replace("_", " "),
                    color="tab:cyan", linewidth=1.3)
        _decorate(ax, title="Average endothelial cell volume", ylabel="Volume (voxels)")
        ax.legend(fontsize=8)

    # growth rate (smoothed)
    if rate_cols:
        ax = axes[panel]; panel += 1
        for col in rate_cols:
            ax.plot(mcs, _smooth(df[col]),
                    label=col.replace("_", " ") + " (smooth)", linewidth=1.4, color="tab:cyan")
        ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.5)
        _decorate(ax, title="Endothelial growth rate (smoothed)",
                  xlabel="MCS", ylabel="ΔVol / MCS")
        ax.legend(fontsize=8)

    fig.tight_layout()
    return fig


# ── 3. Signaling fields ───────────────────────────────────────────────────────

def plot_signaling_fields(df: pd.DataFrame) -> plt.Figure | None:
    """Two panels: oxygen (drop and recovery) · VEGF2 (dual axis with hypoxic fraction)."""
    o2_cols = available_columns(
        df, ["mean_tumor_oxygen"]  # Only tumor oxygen, remove vascular/endothelial
    )
    vegf_cols = available_columns(
        df, ["mean_tumor_vegf2", "mean_tumor_vegf1"]  # Only tumor VEGF, remove vascular
    )
    if not o2_cols and not vegf_cols:
        return None

    nrows = int(bool(o2_cols)) + int(bool(vegf_cols))
    fig, raw_axes = plt.subplots(nrows, 1, figsize=(11, 4.5 * nrows), sharex=True)
    axes = [raw_axes] if nrows == 1 else list(raw_axes)
    fig.suptitle("3 · Signaling Fields", fontsize=12, fontweight="bold")
    mcs = df["mcs"]
    panel = 0

    if o2_cols:
        ax = axes[panel]; panel += 1
        # Only plot tumor oxygen
        if "mean_tumor_oxygen" in df.columns:
            ax.plot(mcs, df["mean_tumor_oxygen"], label="Tumor (mean)", color="tab:red", linestyle="-", linewidth=1.4)
        _decorate(ax,
                  title="Oxygen field\n(tumor consumes O\u2082; angiogenesis restores supply)",
                  ylabel="[O\u2082] (a.u.)")
        ax.legend(fontsize=8)

    if vegf_cols:
        ax = axes[panel]; panel += 1
        # Only plot tumor VEGF fields
        VEGF_STYLE: dict[str, tuple[str, str, str]] = {
            "mean_tumor_vegf2":   ("VEGF2 @ tumor",    "tab:orange", "-"),
            "mean_tumor_vegf1":   ("VEGF1 @ tumor",    "tab:purple", ":"),
        }
        for col, (label, color, ls) in VEGF_STYLE.items():
            if col in df.columns:
                ax.plot(mcs, _smooth(df[col]),
                        label=label + " (smooth)", color=color,
                        linestyle=ls, linewidth=1.4)
        _decorate(ax,
                  title="VEGF fields\n(rise with hypoxia; fall as new vessels supply O\u2082)",
                  xlabel="MCS", ylabel="VEGF (a.u.)")
        ax.legend(loc="upper left", fontsize=8)

    axes[-1].set_xlabel("MCS")
    fig.tight_layout()
    return fig


# ── 4. System-level dynamics ──────────────────────────────────────────────────

def plot_system_dynamics(df: pd.DataFrame) -> plt.Figure | None:
    """Two panels: tumor growth rate (smoothed, annotated) · O2 recovery vs active sprouts."""
    rate_col = next(
        (c for c in ["avg_tumor_volume_growth_rate", "total_tumor_volume_growth_rate"]
         if c in df.columns), None
    )
    o2_col = "mean_tumor_oxygen" if "mean_tumor_oxygen" in df.columns else None
    sprout_col = None  # Remove vascular cell overlays
    if not rate_col and not o2_col:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    fig.suptitle("4 · System-level Dynamics", fontsize=12, fontweight="bold")
    mcs = df["mcs"]

    # ── tumor growth rate ─────────────────────────────────────────────────────
    ax = axes[0]
    if rate_col:
        raw = df[rate_col]
        smooth = _smooth(raw)
        ax.plot(mcs, raw, color="tab:blue", alpha=0.20, linewidth=0.7)
        ax.plot(mcs, smooth, color="tab:blue", linewidth=2.0,
                label="Tumor growth rate (smooth)")
        ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.5)
        # mark zero-crossings (growth → stasis)
        zero_cross = mcs[
            (smooth.shift(1).gt(0) & smooth.le(0)) |
            (smooth.shift(1).lt(0) & smooth.ge(0))
        ]
        crossing_candidates = [int(v) for v in zero_cross.tolist()]
        span = int(mcs.max() - mcs.min()) if len(mcs) else 0
        min_gap = max(10, int(ZERO_CROSSING_MIN_GAP_FRACTION * span))
        for cx in _debounce_crossings(crossing_candidates, min_gap):
            ax.axvline(cx, color="tab:red", linestyle=":", alpha=0.6, linewidth=1.0)
            ax.annotate(
                "rate\u2248 0",
                xy=(cx, 0),
                xytext=(cx + max(5, int((mcs.max() - mcs.min()) * 0.03)),
                        float(smooth.abs().max()) * 0.4),
                fontsize=7, color="tab:red",
                arrowprops=dict(arrowstyle="->", color="tab:red", lw=0.8),
            )
    if "avg_tumor_target_growth_rate" in df.columns:
        ax.plot(mcs, _smooth(df["avg_tumor_target_growth_rate"]),
                color="tab:green", linewidth=1.5, linestyle="--",
                label="Target growth rate (smooth)")
    _decorate(ax,
              title="Tumor growth rate\n(fast under O\u2082; slows with hypoxia; may recover after angiogenesis)",
              ylabel="\u0394Vol / MCS")
    ax.legend(fontsize=8)

    # ── O2 recovery vs active sprout count ───────────────────────────────────
    ax = axes[1]
    if o2_col:
        ax.plot(mcs, df[o2_col], color="tab:cyan", linewidth=1.5,
                label="Mean tumor O\u2082")
        ax.set_ylabel("Mean [O\u2082] (a.u.)")
        ax.tick_params(axis="y", colors="tab:cyan")
        ax.spines["left"].set_color("tab:cyan")
    # Removed overlay of active neovascular cells
    _decorate(ax,
              title="Oxygen recovery vs. active sprout count\n(sprouting precedes O\u2082 restoration; observe the phase lag)",
              xlabel="MCS", ylabel="Mean [O\u2082] (a.u.)")
    if o2_col:
        ax.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    return fig


# ── 5. HIF-1a network (only if columns are present) ───────────────────────────

def plot_hif1a_network(df: pd.DataFrame) -> plt.Figure | None:
    """HIF-1a dynamics and VEGF drive, only emitted when the network was enabled."""
    if "mean_tumor_hif1a" not in df.columns:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    fig.suptitle("5 · HIF-1\u03b1 Gene Network", fontsize=12, fontweight="bold")
    mcs = df["mcs"]

    ax = axes[0]
    ax.plot(mcs, df["mean_tumor_hif1a"], color="tab:purple", linewidth=1.5,
            label="Mean HIF-1\u03b1")
    _decorate(ax, title="Intracellular HIF-1\u03b1 (stabilised by low O\u2082, degraded constitutively)",
              ylabel="HIF-1\u03b1 (a.u.)")
    if "mean_tumor_oxygen" in df.columns:
        ax2 = _twin_right(ax, "tab:red")
        ax2.plot(mcs, df["mean_tumor_oxygen"], color="tab:red",
                 linestyle="--", linewidth=1.1, alpha=0.7, label="Mean O\u2082")
        ax2.set_ylabel("Mean [O\u2082]")
        ax2.legend(loc="upper right", fontsize=8)
    ax.legend(loc="upper left", fontsize=8)

    ax = axes[1]
    if "mean_tumor_vegf_drive" in df.columns:
        ax.plot(mcs, df["mean_tumor_vegf_drive"], color="tab:orange", linewidth=1.5,
                label="Mean VEGF drive")
    _decorate(ax,
              title="Intracellular VEGF drive (Hill function of HIF-1\u03b1)\nboosts effective VEGF2 seen by endothelium",
              xlabel="MCS", ylabel="VEGF drive (a.u.)")
    if "mean_tumor_vegf2" in df.columns:
        ax2 = _twin_right(ax, "sienna")
        ax2.plot(mcs, _smooth(df["mean_tumor_vegf2"]), color="sienna",
                 linestyle="--", linewidth=1.0, alpha=0.7, label="Mean VEGF2 field (smooth)")
        ax2.set_ylabel("Mean VEGF2 field")
        ax2.legend(loc="upper right", fontsize=8)
    ax.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    return fig


def generate_plots(df: pd.DataFrame) -> list[tuple[str, plt.Figure]]:
    figures: list[tuple[str, plt.Figure]] = []
    for stem, figure in [
        ("1_tumor_dynamics",    plot_tumor_dynamics(df)),
        ("2_vascular_response", plot_vascular_response(df)),
        ("3_signaling_fields",  plot_signaling_fields(df)),
        ("4_system_dynamics",   plot_system_dynamics(df)),
        ("5_hif1a_network",     plot_hif1a_network(df)),
    ]:
        if figure is not None:
            figures.append((stem, figure))
    return figures


def main() -> int:
    args = parse_args()
    df = load_metrics(args.input)
    output_dir = resolve_output_dir(args.input, args.output_dir)

    figures = generate_plots(df)
    if not figures:
        raise ValueError("No plottable metric columns were found in the CSV.")

    saved_paths = []
    for stem, figure in figures:
        saved_paths.append(save_figure(figure, output_dir, args.prefix, stem, args.format, args.dpi))

    print("Generated plots:")
    for path in saved_paths:
        print(f" - {path}")

    if args.show:
        reloaded_df = load_metrics(args.input)
        for _, figure in generate_plots(reloaded_df):
            figure.show()
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

