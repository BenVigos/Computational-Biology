"""
Comparison figures across different simulation configurations.

Figure A: Tumor growth (tumor_growth_with_hif1a) vs full model (full_trajectory)
    - Panel a: tumor effective radius vs time
    - Panel b: cell counts (normal, hypoxic, necrotic, total) vs time
    - Panel c: hypoxic fraction vs time

Figure B: Endothelial activation without (late_stage_angiogenesis_only)
          and with (late_stage_angiogenesis_with_hif1a) HIF-1α
    - Panel a: HIF-1α and VEGF vs time
    - Panel b: active and inactive neovascular cell counts vs time
    - Panel c: O2 and minimum endothelial cell distance from tumor vs time

Figure C: Combined tumor dynamics (tumor_growth_with_hif1a vs full_trajectory)
    - Panel a: tumor cell counts vs time
    - Panel b: VEGF and HIF-1α vs time
    - Panel c: O2 and minimum vessel distances vs time
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd

# ── Shared helpers (mirror style of plot_multi_run_angiogenesis.py) ────────────
SMOOTH_WINDOW = 15
TIME_SHIFT_MCS = 300
INIT_LABEL = "Initialization period"

COLORS = {
    "full_trajectory":                   "#1f77b4",   # blue
    "tumor_growth_with_hif1a":           "#d62728",   # red
    "late_stage_angiogenesis_only":      "#2ca02c",   # green
    "late_stage_angiogenesis_with_hif1a":"#ff7f0e",   # orange
}

LABELS = {
    "full_trajectory":                   "With angiogenesis",
    "tumor_growth_with_hif1a":           "Without angiogenesis",
    "late_stage_angiogenesis_only":      "Without HIF-1α",
    "late_stage_angiogenesis_with_hif1a":"With HIF-1α",
}

# Metric style mapping keeps visual semantics consistent across all figures.
METRIC_STYLES = {
    "oxygen": ("-", 1.8),
    "vegf": ("-", 1.8),
    "hif": ("--", 1.5),
    "cell_total": ("-", 2.0),
    "cell_normoxic": ("--", 1.3),
    "cell_hypoxic": ("-.", 1.3),
    "cell_necrotic": (":", 1.3),
    "active_neovascular": ("--", 1.4),
    "inactive_neovascular": (":", 1.4),
    "dist_any_vessel": ("--", 1.5),
    "dist_sprout": (":", 1.5),
    "radius_effective": ("-", 2.0),
    "radius_mean": ("--", 1.6),
}


def _smooth(s: pd.Series, w: int = SMOOTH_WINDOW) -> pd.Series:
    return s.rolling(w, center=True, min_periods=1).mean()


def _get_stats(df: pd.DataFrame, col: str, smooth: bool = False):
    """Return (shifted_mcs_index, mean_series, ci_series) for *col* grouped by mcs."""
    if col not in df.columns:
        return None, None, None
    grp = df.groupby("mcs")[col].agg(["mean", "std", "count"])
    mean = grp["mean"]
    if smooth:
        mean = _smooth(mean)
    ci = 1.96 * grp["std"] / np.sqrt(grp["count"])
    ci = ci.fillna(0)
    shifted_mcs = grp.index - TIME_SHIFT_MCS
    return shifted_mcs, mean, ci


def _plot_ci(ax: plt.Axes, df: pd.DataFrame, col: str, label: str, color: str,
             lw: float = 1.6, ls: str = "-", smooth: bool = False,
             alpha_fill: float = 0.15) -> bool:
    mcs, mean, ci = _get_stats(df, col, smooth=smooth)
    if mcs is None:
        return False
    ax.plot(mcs, mean, label=label, color=color, linewidth=lw, linestyle=ls)
    ax.fill_between(mcs, mean - ci, mean + ci, color=color, alpha=alpha_fill)
    return True


def _decorate(
    ax: plt.Axes,
    *,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    show_initialization: bool = True,
) -> None:
    if title:
        ax.set_title(title, fontsize=9)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if show_initialization:
        # Highlight the pre-growth initialization segment [-300, 0].
        ax.axvspan(-TIME_SHIFT_MCS, 0, color="0.9", alpha=0.45, zorder=0)
        ymin, ymax = ax.get_ylim()
        y_text = ymin + 0.3 * (ymax - ymin)
        ax.text(-TIME_SHIFT_MCS + 6, y_text, INIT_LABEL, fontsize=10, color="0.35", rotation=90)
    # Force major ticks to align on 300-MCS boundaries: -300, 0, 300, ...
    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.grid(alpha=0.25, linestyle="--")


def _twin_right(ax: plt.Axes, color: str) -> plt.Axes:
    ax2 = ax.twinx()
    ax2.spines["right"].set_color(color)
    ax2.tick_params(axis="y", colors=color)
    ax2.yaxis.label.set_color(color)
    return ax2


def _load_config(results_root: Path, config_name: str) -> pd.DataFrame | None:
    """Load and concatenate all CSV runs for a given config folder."""
    config_dir = results_root / config_name
    if not config_dir.exists():
        print(f"  [warn] directory not found: {config_dir}")
        return None
    csvs = sorted(config_dir.glob("angiogenesis_metrics_run*.csv"))
    if not csvs:
        print(f"  [warn] no CSV files found in: {config_dir}")
        return None
    dfs = []
    for p in csvs:
        df = pd.read_csv(p)
        df["config"] = config_name
        dfs.append(df)
    combined = pd.concat(dfs, ignore_index=True)
    return combined


# ══════════════════════════════════════════════════════════════════════════════
# Figure A – Tumor growth vs full model
# ══════════════════════════════════════════════════════════════════════════════
def figure_A_tumor_vs_angiogenesis(
    df_tumor: pd.DataFrame,
    df_full: pd.DataFrame,
) -> plt.Figure:
    """
    Compare tumor growth only (tumor_growth_with_hif1a) against the full model
    (full_trajectory).

    Panel a: tumor effective radius over time
    Panel b: cell counts (total, normoxic, hypoxic, necrotic) over time
    Panel c: hypoxic fraction over time
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(
        "Figure A · Tumor Growth - With and Without Angiogenesis",
        fontsize=12, fontweight="bold",
    )

    configs = [
        ("tumor_growth_with_hif1a", df_tumor),
        ("full_trajectory",         df_full),
    ]

    # ── Panel a: effective + mean radial distance ───────────────────────────
    ax = axes[0]
    for cfg, df in configs:
        eff_ls, eff_lw = METRIC_STYLES["radius_effective"]
        mean_ls, mean_lw = METRIC_STYLES["radius_mean"]
        _plot_ci(ax, df, "tumor_effective_radius", f"{LABELS[cfg]} - Effective radius", COLORS[cfg], lw=eff_lw, ls=eff_ls)
        _plot_ci(ax, df, "tumor_mean_radius", f"{LABELS[cfg]} - Mean radial distance", COLORS[cfg], lw=mean_lw, ls=mean_ls)
    _decorate(ax,
              title="(a) Tumor radius metrics",
              xlabel="MCS", ylabel="Radius (lattice units)")

    config_handles = [
        mlines.Line2D([], [], color=COLORS[c], linewidth=2, label=LABELS[c])
        for c, _ in configs
    ]
    radius_handles = [
        mlines.Line2D([], [], color="black", linewidth=METRIC_STYLES["radius_effective"][1], linestyle=METRIC_STYLES["radius_effective"][0], label="Effective radius"),
        mlines.Line2D([], [], color="black", linewidth=METRIC_STYLES["radius_mean"][1], linestyle=METRIC_STYLES["radius_mean"][0], label="Mean radial distance"),
    ]
    legend_colors = ax.legend(handles=config_handles, loc="upper left", fontsize=7, title="Color")
    legend_styles = ax.legend(handles=radius_handles, loc="upper right", fontsize=7, title="Line style")
    ax.add_artist(legend_colors)

    # ── Panel b: cell counts ─────────────────────────────────────────────────
    ax = axes[1]
    cell_specs = [
        ("tumor_like_cells", "Total", METRIC_STYLES["cell_total"]),
        ("normal_cells", "Normoxic", METRIC_STYLES["cell_normoxic"]),
        ("hypoxic_cells", "Hypoxic", METRIC_STYLES["cell_hypoxic"]),
        ("necrotic_cells", "Necrotic", METRIC_STYLES["cell_necrotic"]),
    ]
    for cfg, df in configs:
        for col, cell_label, (ls, lw) in cell_specs:
            _plot_ci(ax, df, col, f"{LABELS[cfg]} - {cell_label}", COLORS[cfg], lw=lw, ls=ls)

    _decorate(ax,
              title="(b) Cell counts over time",
              xlabel="MCS", ylabel="Number of cells")

    config_handles = [
        mlines.Line2D([], [], color=COLORS[c], linewidth=2, label=LABELS[c])
        for c, _ in configs
    ]
    cell_handles = [
        mlines.Line2D([], [], color="black", linewidth=lw, linestyle=ls, label=cell_label)
        for _, cell_label, (ls, lw) in cell_specs
    ]
    legend_colors = ax.legend(handles=config_handles, loc="upper left", fontsize=7, title="Color")
    legend_styles = ax.legend(handles=cell_handles, loc="upper right", fontsize=7, title="Line style")
    ax.add_artist(legend_colors)

    # ── Panel c: hypoxic fraction ────────────────────────────────────────────
    ax = axes[2]
    for cfg, df in configs:
        _plot_ci(
            ax,
            df,
            "hypoxic_fraction",
            LABELS[cfg],
            COLORS[cfg],
            lw=1.9,
            ls=METRIC_STYLES["cell_hypoxic"][0],
            smooth=True,
        )
    _decorate(
        ax,
        title="(c) Hypoxic fraction over time",
        xlabel="MCS",
        ylabel="Hypoxic fraction",
    )
    ax.legend(fontsize=8)

    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Figure B – Endothelial activation ± HIF-1α
# ══════════════════════════════════════════════════════════════════════════════
def figure_B_endothelial_activation(
    df_no_hif: pd.DataFrame,
    df_hif: pd.DataFrame,
) -> plt.Figure:
    """
    Compare late-stage angiogenesis without HIF-1a (late_stage_angiogenesis_only)
    vs with HIF-1a (late_stage_angiogenesis_with_hif1a).

    Panel a: HIF-1a and VEGF over time
    Panel b: active and inactive neovascular cell counts over time
    Panel c: O2 and minimum endothelial-cell distance from tumor over time
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)
    fig.suptitle(
        "Figure B · Endothelial Activation - With and Without HIF-1α",
        fontsize=12, fontweight="bold",
    )

    configs = [
        ("late_stage_angiogenesis_only",       df_no_hif),
        ("late_stage_angiogenesis_with_hif1a", df_hif),
    ]

    # ── Panel a: HIF-1a and VEGF ────────────────────────────────────────────
    ax = axes[0]
    for cfg, df in configs:
        vegf_ls, vegf_lw = METRIC_STYLES["vegf"]
        _plot_ci(ax, df, "mean_tumor_vegf2", f"VEGF - {LABELS[cfg]}", COLORS[cfg], lw=vegf_lw, ls=vegf_ls, smooth=True)

    ax.set_ylabel("Mean VEGF (a.u.)")
    _decorate(ax, title="(a) VEGF and HIF-1a over time")
    ax.legend(loc="upper left", fontsize=8)

    # HIF-1a on right axis (distinct by config color + style)
    ax2 = _twin_right(ax, "0.25")
    hif_handles = []
    for cfg, df in configs:
        hif_ls, hif_lw = METRIC_STYLES["hif"]
        hif_ls_cfg = hif_ls if cfg == "late_stage_angiogenesis_with_hif1a" else ":"
        plotted = _plot_ci(ax2, df, "mean_tumor_hif1a", f"HIF-1a - {LABELS[cfg]}", COLORS[cfg], lw=hif_lw, ls=hif_ls_cfg)
        if plotted:
            hif_handles.append(
                mlines.Line2D([], [], color=COLORS[cfg], linewidth=hif_lw, linestyle=hif_ls_cfg, label=f"HIF-1a - {LABELS[cfg]}")
            )
    ax2.set_ylabel("Mean HIF-1a (a.u.)")
    if hif_handles:
        ax2.legend(handles=hif_handles, loc="upper right", fontsize=8)

    # ── Panel b: active/inactive neovascular only ───────────────────────────
    ax = axes[1]
    for cfg, df in configs:
        act_ls, act_lw = METRIC_STYLES["active_neovascular"]
        inact_ls, inact_lw = METRIC_STYLES["inactive_neovascular"]
        _plot_ci(ax, df, "active_neovascular_cells", f"Active neovascular - {LABELS[cfg]}", COLORS[cfg], lw=act_lw, ls=act_ls)
        _plot_ci(ax, df, "inactive_neovascular_cells", f"Inactive neovascular - {LABELS[cfg]}", COLORS[cfg], lw=inact_lw, ls=inact_ls)

    _decorate(ax, title="(b) Active vs inactive neovascular cells", ylabel="Cell count")
    config_handles = [
        mlines.Line2D([], [], color=COLORS[c], linewidth=2, label=LABELS[c])
        for c, _ in configs
    ]
    type_handles = [
        mlines.Line2D([], [], color="black", linewidth=METRIC_STYLES["active_neovascular"][1], linestyle=METRIC_STYLES["active_neovascular"][0], label="Active neovascular"),
        mlines.Line2D([], [], color="black", linewidth=METRIC_STYLES["inactive_neovascular"][1], linestyle=METRIC_STYLES["inactive_neovascular"][0], label="Inactive neovascular"),
    ]
    legend_colors = ax.legend(handles=config_handles, loc="upper left", fontsize=8, title="Color")
    legend_styles = ax.legend(handles=type_handles, loc="upper right", fontsize=8, title="Line style")
    ax.add_artist(legend_colors)

    # ── Panel c: O2 + minimum distance to vessel ────────────────────────────
    ax = axes[2]
    for cfg, df in configs:
        o2_ls, o2_lw = METRIC_STYLES["oxygen"]
        _plot_ci(ax, df, "mean_tumor_oxygen", f"O2 - {LABELS[cfg]}", COLORS[cfg], lw=o2_lw, ls=o2_ls)

    ax.set_ylabel("Mean O2 (a.u.)")
    _decorate(ax, title="(c) O2 and minimum vessel distance from tumor", xlabel="MCS")
    ax.legend(loc="upper left", fontsize=8)

    ax2 = _twin_right(ax, "0.25")
    for cfg, df in configs:
        dv_ls, dv_lw = METRIC_STYLES["dist_any_vessel"]
        ds_ls, ds_lw = METRIC_STYLES["dist_sprout"]
        # _plot_ci(ax2, df, "min_tumor_to_vessel_distance", f"Min dist to any vessel - {LABELS[cfg]}", COLORS[cfg], lw=dv_lw, ls=dv_ls)
        _plot_ci(ax2, df, "min_dist_to_sprout", f"Min dist to sprout - {LABELS[cfg]}", COLORS[cfg], lw=ds_lw, ls=ds_ls)
    ax2.set_ylabel("Distance (lattice units)")
    ax2.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Figure C – Combined tumor dynamics
# ══════════════════════════════════════════════════════════════════════════════
def figure_C_combined_tumor_dynamics(
    df_tumor: pd.DataFrame,
    df_full: pd.DataFrame,
) -> plt.Figure:
    """
    Compare tumor_growth_with_hif1a vs full_trajectory across:

    Panel a: tumor cell counts (total, normoxic, hypoxic, necrotic) over time
    Panel b: VEGF and HIF-1a over time
    Panel c: O2 and minimum vessel distances over time
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)
    fig.suptitle(
        "Figure C · Combined Tumor Dynamics  –  With and Without Angiogenesis",
        fontsize=12, fontweight="bold",
    )

    configs = [
        ("tumor_growth_with_hif1a", df_tumor),
        ("full_trajectory",         df_full),
    ]
    cell_specs = [
        ("tumor_like_cells", "Total", METRIC_STYLES["cell_total"]),
        ("normal_cells", "Normoxic", METRIC_STYLES["cell_normoxic"]),
        ("hypoxic_cells", "Hypoxic", METRIC_STYLES["cell_hypoxic"]),
        ("necrotic_cells", "Necrotic", METRIC_STYLES["cell_necrotic"]),
    ]

    # ── Panel a: cell counts ─────────────────────────────────────────────────
    ax = axes[0]
    for cfg, df in configs:
        for col, cell_label, (ls, lw) in cell_specs:
            _plot_ci(ax, df, col, f"{LABELS[cfg]} - {cell_label}", COLORS[cfg], lw=lw, ls=ls)
    _decorate(ax, title="(a) Tumor cell counts over time", ylabel="Number of cells")

    config_handles = [
        mlines.Line2D([], [], color=COLORS[c], linewidth=2, label=LABELS[c])
        for c, _ in configs
    ]
    cell_handles = [
        mlines.Line2D([], [], color="black", linewidth=lw, linestyle=ls, label=cell_label)
        for _, cell_label, (ls, lw) in cell_specs
    ]
    legend_colors = ax.legend(handles=config_handles, loc="upper left", fontsize=7, title="Color")
    legend_styles = ax.legend(handles=cell_handles, loc="upper right", fontsize=7, title="Line style")
    ax.add_artist(legend_colors)

    # ── Panel b: VEGF and HIF-1a ─────────────────────────────────────────────
    ax = axes[1]
    for cfg, df in configs:
        vegf_ls, vegf_lw = METRIC_STYLES["vegf"]
        _plot_ci(ax, df, "mean_tumor_vegf2", f"VEGF - {LABELS[cfg]}", COLORS[cfg], lw=vegf_lw, ls=vegf_ls, smooth=True)
    ax.set_ylabel("Mean VEGF (a.u.)")
    _decorate(ax, title="(b) VEGF and HIF-1a over time")
    ax.legend(loc="upper left", fontsize=8)

    ax2 = _twin_right(ax, "0.25")
    for cfg, df in configs:
        hif_ls, hif_lw = METRIC_STYLES["hif"]
        hif_ls_cfg = hif_ls if cfg == "full_trajectory" else ":"
        _plot_ci(ax2, df, "mean_tumor_hif1a", f"HIF-1a - {LABELS[cfg]}", COLORS[cfg], lw=hif_lw, ls=hif_ls_cfg)
    ax2.set_ylabel("Mean HIF-1a (a.u.)")
    ax2.legend(loc="upper right", fontsize=8)

    # ── Panel c: O2 + minimum distances ─────────────────────────────────────
    ax = axes[2]
    for cfg, df in configs:
        o2_ls, o2_lw = METRIC_STYLES["oxygen"]
        _plot_ci(ax, df, "mean_tumor_oxygen", f"O2 - {LABELS[cfg]}", COLORS[cfg], lw=o2_lw, ls=o2_ls)
    _decorate(
        ax,
        title="(c) Tumor oxygen and minimum vessel distance",
        xlabel="MCS",
        ylabel="Mean O2 (a.u.)",
    )
    ax.legend(loc="upper left", fontsize=8)

    ax2 = _twin_right(ax, "0.25")
    for cfg, df in configs:
        dv_ls, dv_lw = METRIC_STYLES["dist_any_vessel"]
        ds_ls, ds_lw = METRIC_STYLES["dist_sprout"]
        _plot_ci(ax2, df, "min_tumor_to_vessel_distance", f"Min dist to any vessel - {LABELS[cfg]}", COLORS[cfg], lw=dv_lw, ls=dv_ls)
        # _plot_ci(ax2, df, "min_dist_to_sprout", f"Min dist to neovascular sprout - {LABELS[cfg]}", COLORS[cfg], lw=ds_lw, ls=ds_ls)
    ax2.set_ylabel("Distance (lattice units)")
    ax2.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# CLI entry point
# ══════════════════════════════════════════════════════════════════════════════
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate comparison figures across simulation configurations."
    )
    p.add_argument(
        "--results-dir", type=Path,
        default=Path(__file__).parent.parent / "results",
        help="Root folder containing per-config sub-directories (default: ../results)",
    )
    p.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).parent.parent / "results" / "comparison_figures",
        help="Where to save the output PNGs",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    root = args.results_dir

    print("Loading data …")
    df_full   = _load_config(root, "full_trajectory")
    df_tumor  = _load_config(root, "tumor_growth_with_hif1a")
    df_no_hif = _load_config(root, "late_stage_angiogenesis_only")
    df_hif    = _load_config(root, "late_stage_angiogenesis_with_hif1a")

    figures = []

    if df_tumor is not None and df_full is not None:
        figures.append(("figure_A_tumor_vs_angiogenesis",
                         figure_A_tumor_vs_angiogenesis(df_tumor, df_full)))
        figures.append(("figure_C_combined_tumor_dynamics",
                         figure_C_combined_tumor_dynamics(df_tumor, df_full)))
    else:
        print("[warn] Skipping figures A and C – missing data.")

    if df_no_hif is not None and df_hif is not None:
        figures.append(("figure_B_endothelial_activation",
                         figure_B_endothelial_activation(df_no_hif, df_hif)))
    else:
        print("[warn] Skipping figure B – missing data.")

    for stem, fig in figures:
        out = args.output_dir / f"{stem}.png"
        fig.savefig(out, dpi=160, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved → {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

