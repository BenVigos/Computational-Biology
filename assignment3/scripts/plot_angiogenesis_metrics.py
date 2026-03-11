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


DEFAULT_INPUT = Path("../Angiogenesis/results/angiogenesis_metrics.csv")
DEFAULT_FORMAT = "png"
DEFAULT_DPI = 160


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


def plot_counts_and_fraction(df: pd.DataFrame) -> plt.Figure | None:
    count_columns = available_columns(
        df,
        [
            "tumor_cells",
            "normal_cells",
            "hypoxic_cells",
            "necrotic_cells",
            "tumor_like_cells",
            "endothelial_cells",
            "vascular_like_cells",
            "active_neovascular_cells",
            "vascular_cells",
            "inactive_neovascular_cells",
        ],
    )
    fraction_columns = available_columns(df, ["hypoxic_fraction", "necrotic_fraction"])
    if not count_columns and not fraction_columns:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    if count_columns:
        for column in count_columns:
            axes[0].plot(df["mcs"], df[column], label=column.replace("_", " "))
        axes[0].set_ylabel("Cell count")
        axes[0].set_title("Cell populations over time")
        axes[0].grid(alpha=0.3)
        axes[0].legend()
    else:
        axes[0].set_visible(False)

    if fraction_columns:
        for column in fraction_columns:
            axes[1].plot(df["mcs"], df[column], label=column.replace("_", " "))
        axes[1].set_ylabel("Fraction")
        axes[1].set_xlabel("MCS")
        axes[1].set_title("Tumor phenotype fractions over time")
        axes[1].grid(alpha=0.3)
        axes[1].legend()
    else:
        axes[1].set_visible(False)

    return fig


def plot_tumor_volume(df: pd.DataFrame) -> plt.Figure | None:
    columns = available_columns(
        df,
        [
            "avg_tumor_volume",
            "total_tumor_volume",
            "avg_tumor_target_volume",
            "avg_vascular_volume",
            "avg_endothelial_volume",
            "mean_hif",
        ],
    )
    if not columns:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    upper_columns = available_columns(
        df,
        [
            "avg_tumor_volume",
            "avg_tumor_target_volume",
            "total_tumor_volume",
            "avg_vascular_volume",
            "avg_endothelial_volume",
        ],
    )
    for column in upper_columns:
        axes[0].plot(df["mcs"], df[column], label=column.replace("_", " "))
    axes[0].set_ylabel("Volume")
    axes[0].set_title("Tumor size metrics")
    axes[0].grid(alpha=0.3)
    if upper_columns:
        axes[0].legend()

    if "mean_hif" in df.columns:
        axes[1].plot(df["mcs"], df["mean_hif"], color="tab:purple", label="mean HIF")
        axes[1].set_ylabel("Mean HIF")
        axes[1].set_title("Hypoxia signalling")
        axes[1].grid(alpha=0.3)
        axes[1].legend()
    else:
        axes[1].set_visible(False)

    axes[1].set_xlabel("MCS")
    return fig


def plot_growth_rates(df: pd.DataFrame) -> plt.Figure | None:
    columns = available_columns(
        df,
        [
            "avg_tumor_volume_growth_rate",
            "avg_tumor_target_growth_rate",
            "total_tumor_volume_growth_rate",
            "avg_endothelial_volume_growth_rate",
            "avg_vascular_volume_growth_rate",
        ],
    )
    if not columns:
        return None

    fig, ax = plt.subplots(figsize=(10, 5.5))
    for column in columns:
        ax.plot(df["mcs"], df[column], label=column.replace("_", " "))
    ax.axhline(0.0, color="black", linewidth=1, alpha=0.5)
    ax.set_xlabel("MCS")
    ax.set_ylabel("Rate per MCS")
    ax.set_title("Growth-rate metrics")
    ax.grid(alpha=0.3)
    ax.legend()
    return fig


def plot_field_means(df: pd.DataFrame) -> plt.Figure | None:
    tumor_field_columns = available_columns(df, ["mean_tumor_oxygen", "mean_tumor_vegf", "mean_tumor_vegf1", "mean_tumor_vegf2"])
    endothelial_field_columns = available_columns(
        df,
        [
            "mean_endothelial_oxygen",
            "mean_endothelial_vegf",
            "mean_vascular_oxygen",
            "mean_vascular_vegf1",
            "mean_vascular_vegf2",
        ],
    )
    if not tumor_field_columns and not endothelial_field_columns:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    if tumor_field_columns:
        for column in tumor_field_columns:
            axes[0].plot(df["mcs"], df[column], label=column.replace("_", " "))
        axes[0].set_ylabel("Field value")
        axes[0].set_title("Tumor-local field means")
        axes[0].grid(alpha=0.3)
        axes[0].legend()
    else:
        axes[0].set_visible(False)

    if endothelial_field_columns:
        for column in endothelial_field_columns:
            axes[1].plot(df["mcs"], df[column], label=column.replace("_", " "))
        axes[1].set_xlabel("MCS")
        axes[1].set_ylabel("Field value")
        axes[1].set_title("Endothelial-local field means")
        axes[1].grid(alpha=0.3)
        axes[1].legend()
    else:
        axes[1].set_visible(False)

    return fig


def plot_endothelial_metrics(df: pd.DataFrame) -> plt.Figure | None:
    columns = available_columns(
        df,
        [
            "endothelial_cells",
            "vascular_like_cells",
            "active_neovascular_cells",
            "vascular_cells",
            "inactive_neovascular_cells",
            "avg_endothelial_volume",
            "avg_vascular_volume",
        ],
    )
    if not columns:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    population_columns = available_columns(
        df,
        ["endothelial_cells", "vascular_like_cells", "active_neovascular_cells", "vascular_cells", "inactive_neovascular_cells"],
    )
    if population_columns:
        for column in population_columns:
            axes[0].plot(df["mcs"], df[column], label=column.replace("_", " "))
        axes[0].set_ylabel("Cell count")
        axes[0].set_title("Vascular / neovascular populations")
        axes[0].grid(alpha=0.3)
        axes[0].legend()
    else:
        axes[0].set_visible(False)

    size_columns = available_columns(df, ["avg_endothelial_volume", "avg_vascular_volume"])
    if size_columns:
        for column in size_columns:
            axes[1].plot(df["mcs"], df[column], label=column.replace("_", " "))
        axes[1].set_xlabel("MCS")
        axes[1].set_ylabel("Volume")
        axes[1].set_title("Vascular / neovascular size")
        axes[1].grid(alpha=0.3)
        axes[1].legend()
    else:
        axes[1].set_visible(False)

    return fig


def generate_plots(df: pd.DataFrame) -> list[tuple[str, plt.Figure]]:
    figures: list[tuple[str, plt.Figure]] = []
    for stem, figure in [
        ("counts", plot_counts_and_fraction(df)),
        ("tumor_volume", plot_tumor_volume(df)),
        ("growth_rates", plot_growth_rates(df)),
        ("field_means", plot_field_means(df)),
        ("endothelial", plot_endothelial_metrics(df)),
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

