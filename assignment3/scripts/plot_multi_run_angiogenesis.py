"""Create grouped summary plots with 95% CI from multiple angiogenesis monitoring CSVs."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

SMOOTH_WINDOW = 15  
ZERO_CROSSING_MIN_GAP_FRACTION = 0.07

def _smooth(series: pd.Series, window: int = SMOOTH_WINDOW) -> pd.Series:
    return series.rolling(window, center=True, min_periods=1).mean()

def _twin_right(ax: plt.Axes, color: str) -> plt.Axes:
    ax2 = ax.twinx()
    ax2.spines["right"].set_color(color)
    ax2.tick_params(axis="y", colors=color)
    ax2.yaxis.label.set_color(color)
    return ax2

def _decorate(ax: plt.Axes, *, title: str = "", xlabel: str = "", ylabel: str = "") -> None:
    if title: ax.set_title(title, fontsize=9)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25, linestyle="--")

def _debounce_crossings(crossings: list[int], min_gap: int) -> list[int]:
    if not crossings: return []
    min_gap = max(1, int(min_gap))
    kept = [int(crossings[0])]
    for cx in crossings[1:]:
        if int(cx) - kept[-1] >= min_gap:
            kept.append(int(cx))
    return kept

def get_stats(df: pd.DataFrame, col: str, smooth: bool = False):
    """Calculates mean and 95% CI for a given column grouped by mcs."""
    if col not in df.columns:
        return None, None, None
    grouped = df.groupby('mcs')[col].agg(['mean', 'std', 'count'])
    mcs = grouped.index
    mean = grouped['mean']
    if smooth:
        mean = _smooth(mean)
    ci = 1.96 * grouped['std'] / np.sqrt(grouped['count'])
    ci = ci.fillna(0)
    return mcs, mean, ci

def plot_with_ci(ax, df, col, label, color, lw=1.5, ls='-', smooth=False, twin_ax=None):
    """Plots the mean line and shaded 95% CI."""
    mcs, mean, ci = get_stats(df, col, smooth)
    if mcs is None: return False
    
    target_ax = twin_ax if twin_ax else ax
    target_ax.plot(mcs, mean, label=label, color=color, linewidth=lw, linestyle=ls)
    target_ax.fill_between(mcs, mean - ci, mean + ci, color=color, alpha=0.15)
    return True

def available_columns(df: pd.DataFrame, columns: Iterable[str]) -> list[str]:
    return [column for column in columns if column in df.columns]

# ── 1. Tumor dynamics ─────────────────────────────────────────────────────────
def plot_tumor_dynamics(df: pd.DataFrame, config_name: str) -> plt.Figure | None:
    has_counts = any(c in df.columns for c in ["tumor_like_cells", "normal_cells", "hypoxic_cells", "necrotic_cells"])
    has_fractions = any(c in df.columns for c in ["hypoxic_fraction", "necrotic_fraction"])
    has_volumes = any(c in df.columns for c in ["total_tumor_volume", "avg_tumor_volume"])
    if not any([has_counts, has_fractions, has_volumes]): return None

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    fig.suptitle(f"[{config_name}] 1 · Tumor Dynamics", fontsize=12, fontweight="bold")

    # counts
    ax = axes[0, 0]
    plotted = False
    for col, label, color, lw in [
        ("tumor_like_cells", "Total tumor", "black", 2.0),
        ("normal_cells", "Normal", "tab:blue", 1.2),
        ("hypoxic_cells", "Hypoxic", "tab:orange", 1.2),
        ("necrotic_cells", "Necrotic", "tab:red", 1.2),
    ]:
        if plot_with_ci(ax, df, col, label, color, lw): plotted = True
    _decorate(ax, title="Cell counts\n(slow growth → acceleration → O₂ limitation)", ylabel="Number of cells")
    if plotted: ax.legend(fontsize=8)

    # fractions
    ax = axes[0, 1]
    frac_plotted = False
    for col, label, color in [("hypoxic_fraction", "Hypoxic fraction", "tab:orange"), ("necrotic_fraction", "Necrotic fraction", "tab:red")]:
        if plot_with_ci(ax, df, col, label, color, 1.4): frac_plotted = True
    _decorate(ax, title="Phenotype fractions & VEGF2\n(hypoxia drives VEGF2 → angiogenesis → O₂ recovery)", ylabel="Fraction")
    if frac_plotted: ax.set_ylim(bottom=0); ax.legend(loc="upper left", fontsize=8)
    
    if "mean_tumor_vegf2" in df.columns:
        ax2 = _twin_right(ax, "tab:purple")
        plot_with_ci(ax, df, "mean_tumor_vegf2", "Mean VEGF2 (smooth)", "tab:purple", 1.2, '--', smooth=True, twin_ax=ax2)
        ax2.set_ylabel("Mean VEGF2")
        ax2.legend(loc="upper right", fontsize=8)

    # volumes
    ax = axes[1, 0]
    plot_with_ci(ax, df, "total_tumor_volume", "Total volume", "tab:blue", 1.5)
    ax.set_ylabel("Total volume (voxels)")
    if "avg_tumor_volume" in df.columns:
        ax2 = _twin_right(ax, "tab:cyan")
        plot_with_ci(ax, df, "avg_tumor_volume", "Avg cell volume", "tab:cyan", 1.2, '--', twin_ax=ax2)
        ax2.set_ylabel("Avg cell volume (voxels)")
        ax2.legend(loc="upper left", fontsize=8)

    # pressure
    ax = axes[1, 1]
    mcs, target_mean, _ = get_stats(df, "avg_tumor_target_volume")
    _, actual_mean, _ = get_stats(df, "avg_tumor_volume")
    
    if mcs is not None and actual_mean is not None:
        ax.plot(mcs, target_mean, color="tab:green", linewidth=1.5, label="Avg target volume")
        ax.plot(mcs, actual_mean, color="tab:blue", linestyle="--", linewidth=1.2, label="Avg actual volume")
        ax.fill_between(mcs, actual_mean, target_mean, alpha=0.15, color="tab:green", label="Growth pressure")
        ax.legend(loc="upper left", fontsize=8)

    _decorate(ax, title="Growth pressure (target − actual)\n+ HIF-1α if network enabled", xlabel="MCS", ylabel="Volume (voxels)")
    
    if "mean_tumor_hif1a" in df.columns:
        ax2 = _twin_right(ax, "tab:purple")
        plot_with_ci(ax, df, "mean_tumor_hif1a", "Mean HIF-1α", "tab:purple", 1.2, ':', twin_ax=ax2)
        ax2.set_ylabel("Mean HIF-1α")
        ax2.legend(loc="lower right", fontsize=8)

    fig.tight_layout()
    return fig

# ── 2. Vascular response ──────────────────────────────────────────────────────
def plot_vascular_response(df: pd.DataFrame, config_name: str) -> plt.Figure | None:
    pop_cols = available_columns(df, ["active_neovascular_cells", "inactive_neovascular_cells"])
    size_cols = available_columns(df, ["avg_endothelial_volume"])
    rate_cols = available_columns(df, ["avg_endothelial_volume_growth_rate"])
    if not pop_cols and not size_cols and not rate_cols: return None

    n_panels = sum([bool(pop_cols), bool(size_cols), bool(rate_cols)])
    fig, axes = plt.subplots(n_panels, 1, figsize=(11, 3.5 * n_panels), sharex=True)
    axes = [axes] if n_panels == 1 else list(axes)
    fig.suptitle(f"[{config_name}] 2 · Endothelial Response", fontsize=12, fontweight="bold")
    panel = 0

    if pop_cols:
        ax = axes[panel]; panel += 1
        for col, label, color in [("active_neovascular_cells", "Active (tip cell)", "tab:green"), ("inactive_neovascular_cells", "Inactive (stalk)", "tab:olive")]:
            plot_with_ci(ax, df, col, label, color, 1.3)
        _decorate(ax, title="Endothelial populations\n(VEGF2 chemotaxis activates sprouting)", ylabel="Cell count")
        ax.legend(fontsize=8)

    if size_cols:
        ax = axes[panel]; panel += 1
        for col in size_cols: plot_with_ci(ax, df, col, col.replace("_", " "), "tab:cyan", 1.3)
        _decorate(ax, title="Average endothelial cell volume", ylabel="Volume (voxels)")
        ax.legend(fontsize=8)

    if rate_cols:
        ax = axes[panel]; panel += 1
        for col in rate_cols: plot_with_ci(ax, df, col, col.replace("_", " ") + " (smooth)", "tab:cyan", 1.4, smooth=True)
        ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.5)
        _decorate(ax, title="Endothelial growth rate (smoothed)", xlabel="MCS", ylabel="ΔVol / MCS")
        ax.legend(fontsize=8)

    fig.tight_layout()
    return fig

# ── 3. Signaling fields ───────────────────────────────────────────────────────
def plot_signaling_fields(df: pd.DataFrame, config_name: str) -> plt.Figure | None:
    o2_cols = available_columns(df, ["mean_tumor_oxygen"])
    vegf_cols = available_columns(df, ["mean_tumor_vegf2", "mean_tumor_vegf1"])
    if not o2_cols and not vegf_cols: return None

    nrows = int(bool(o2_cols)) + int(bool(vegf_cols))
    fig, axes = plt.subplots(nrows, 1, figsize=(11, 4.5 * nrows), sharex=True)
    axes = [axes] if nrows == 1 else list(axes)
    fig.suptitle(f"[{config_name}] 3 · Signaling Fields", fontsize=12, fontweight="bold")
    panel = 0

    if o2_cols:
        ax = axes[panel]; panel += 1
        plot_with_ci(ax, df, "mean_tumor_oxygen", "Tumor (mean)", "tab:red", 1.4)
        _decorate(ax, title="Oxygen field\n(tumor consumes O₂; angiogenesis restores supply)", ylabel="[O₂] (a.u.)")
        ax.legend(fontsize=8)

    if vegf_cols:
        ax = axes[panel]; panel += 1
        plot_with_ci(ax, df, "mean_tumor_vegf2", "VEGF2 @ tumor (smooth)", "tab:orange", 1.4, smooth=True)
        plot_with_ci(ax, df, "mean_tumor_vegf1", "VEGF1 @ tumor (smooth)", "tab:purple", 1.4, ':', smooth=True)
        _decorate(ax, title="VEGF fields\n(rise with hypoxia; fall as new vessels supply O₂)", xlabel="MCS", ylabel="VEGF (a.u.)")
        ax.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    return fig

# ── 4. System-level dynamics ──────────────────────────────────────────────────
def plot_system_dynamics(df: pd.DataFrame, config_name: str) -> plt.Figure | None:
    rate_col = next((c for c in ["avg_tumor_volume_growth_rate", "total_tumor_volume_growth_rate"] if c in df.columns), None)
    o2_col = "mean_tumor_oxygen" if "mean_tumor_oxygen" in df.columns else None
    if not rate_col and not o2_col: return None

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    fig.suptitle(f"[{config_name}] 4 · System-level Dynamics", fontsize=12, fontweight="bold")

    ax = axes[0]
    if rate_col:
        mcs, mean_rate, _ = get_stats(df, rate_col, smooth=True)
        plot_with_ci(ax, df, rate_col, "Tumor growth rate (smooth)", "tab:blue", 2.0, smooth=True)
        ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.5)
        
        # Debounce zero crossings using the mean
        if mean_rate is not None:
            zero_cross = mcs[(mean_rate.shift(1).gt(0) & mean_rate.le(0)) | (mean_rate.shift(1).lt(0) & mean_rate.ge(0))]
            span = int(mcs.max() - mcs.min()) if len(mcs) else 0
            for cx in _debounce_crossings([int(v) for v in zero_cross.tolist()], max(10, int(ZERO_CROSSING_MIN_GAP_FRACTION * span))):
                ax.axvline(cx, color="tab:red", linestyle=":", alpha=0.6, linewidth=1.0)
                ax.annotate("rate≈ 0", xy=(cx, 0), xytext=(cx + max(5, int(span * 0.03)), float(mean_rate.abs().max()) * 0.4),
                            fontsize=7, color="tab:red", arrowprops=dict(arrowstyle="->", color="tab:red", lw=0.8))

    if "avg_tumor_target_growth_rate" in df.columns:
        plot_with_ci(ax, df, "avg_tumor_target_growth_rate", "Target growth rate (smooth)", "tab:green", 1.5, '--', smooth=True)
    
    _decorate(ax, title="Tumor growth rate\n(fast under O₂; slows with hypoxia; may recover after angiogenesis)", ylabel="ΔVol / MCS")
    ax.legend(fontsize=8)

    ax = axes[1]
    if o2_col:
        plot_with_ci(ax, df, o2_col, "Mean tumor O₂", "tab:cyan", 1.5)
        ax.set_ylabel("Mean [O₂] (a.u.)")
        ax.tick_params(axis="y", colors="tab:cyan")
        ax.spines["left"].set_color("tab:cyan")
    
    _decorate(ax, title="Oxygen recovery\n(sprouting precedes O₂ restoration; observe the phase lag)", xlabel="MCS", ylabel="Mean [O₂] (a.u.)")
    if o2_col: ax.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    return fig

# ── 5. HIF-1a network ─────────────────────────────────────────────────────────
def plot_hif1a_network(df: pd.DataFrame, config_name: str) -> plt.Figure | None:
    if "mean_tumor_hif1a" not in df.columns: return None

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    fig.suptitle(f"[{config_name}] 5 · HIF-1α Gene Network", fontsize=12, fontweight="bold")

    ax = axes[0]
    plot_with_ci(ax, df, "mean_tumor_hif1a", "Mean HIF-1α", "tab:purple", 1.5)
    _decorate(ax, title="Intracellular HIF-1α (stabilised by low O₂, degraded constitutively)", ylabel="HIF-1α (a.u.)")
    
    if "mean_tumor_oxygen" in df.columns:
        ax2 = _twin_right(ax, "tab:red")
        plot_with_ci(ax, df, "mean_tumor_oxygen", "Mean O₂", "tab:red", 1.1, '--', twin_ax=ax2)
        ax2.set_ylabel("Mean [O₂]")
        ax2.legend(loc="upper right", fontsize=8)
    ax.legend(loc="upper left", fontsize=8)

    ax = axes[1]
    if "mean_tumor_vegf_drive" in df.columns:
        plot_with_ci(ax, df, "mean_tumor_vegf_drive", "Mean VEGF drive", "tab:orange", 1.5)
    _decorate(ax, title="Intracellular VEGF drive (Hill function of HIF-1α)\nboosts effective VEGF2 seen by endothelium", xlabel="MCS", ylabel="VEGF drive (a.u.)")
    
    if "mean_tumor_vegf2" in df.columns:
        ax2 = _twin_right(ax, "sienna")
        plot_with_ci(ax, df, "mean_tumor_vegf2", "Mean VEGF2 field (smooth)", "sienna", 1.0, '--', smooth=True, twin_ax=ax2)
        ax2.set_ylabel("Mean VEGF2 field")
        ax2.legend(loc="upper right", fontsize=8)
    ax.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()

def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load and combine all CSVs
    dfs = []
    for path in args.inputs:
        if path.exists():
            df = pd.read_csv(path)
            df["config"] = path.parent.name
            dfs.append(df)
            
    if not dfs:
        raise ValueError("No valid CSV files found.")
    
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Generate the 5 dashboard plots for EACH configuration
    for config_name, group_df in combined_df.groupby("config"):
        print(f"Generating dashboard for: {config_name}")
        
        # Create a specific directory for this config to store plots alongside CSVs
        config_out_dir = args.output_dir / config_name
        config_out_dir.mkdir(parents=True, exist_ok=True)
        
        figures = [
            ("1_tumor_dynamics", plot_tumor_dynamics(group_df, config_name)),
            ("2_vascular_response", plot_vascular_response(group_df, config_name)),
            ("3_signaling_fields", plot_signaling_fields(group_df, config_name)),
            ("4_system_dynamics", plot_system_dynamics(group_df, config_name)),
            ("5_hif1a_network", plot_hif1a_network(group_df, config_name)),
        ]
        
        for stem, fig in figures:
            if fig is not None:
                out_path = config_out_dir / f"{stem}.png"
                fig.savefig(out_path, dpi=160, bbox_inches="tight")
                plt.close(fig)
                print(f" - Saved {out_path.name}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())