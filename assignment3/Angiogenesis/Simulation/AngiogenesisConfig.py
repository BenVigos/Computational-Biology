from dataclasses import dataclass, replace


@dataclass(frozen=True)
class ModelConfig:
    preset_name: str = "paper_full"

    # Runtime toggles
    enable_initial_vascular_strip: bool = True
    enable_initial_sprouts: bool = True
    enable_initial_tumor_mass: bool = True
    enable_type_switching: bool = True
    enable_tumor_growth: bool = True
    enable_vascular_growth: bool = True
    enable_mitosis: bool = True
    enable_python_monitoring: bool = True

    # Timescale separation: run diffusion for N MCS, then growth/mitosis for G MCS.
    enable_timescale_separation: bool = True
    diffusion_relaxation_mcs: int = 100
    growth_window_mcs: int = 10
    timescale_cycle_offset_mcs: int = 0

    # 2D lattice / runtime controls
    lattice_x: int = 100
    lattice_y: int = 100
    steps: int = 2000

    # Relative placement controls (fractions of lattice size)
    vessel_boundary_thickness_fraction: float = 0.05
    vessel_boundary_margin_fraction: float = 0.0
    tip_cell_x_offset_fraction: float = 0.0
    tip_cell_width_fraction: float = 0.06
    tip_cell_height_fraction: float = 0.06
    tip_cell_center_y_fractions: tuple[float, ...] = (0.22, 0.48, 0.74)
    tumor_x_min_fraction: float = 0.48
    tumor_x_max_fraction: float = 0.52
    tumor_y_min_fraction: float = 0.49
    tumor_y_max_fraction: float = 0.51
    tumor_seed_size_fraction: float = 0.1

    # Paper-derived phenotype thresholds and timing
    area_thresh: float = 1.0
    nutrient_thresh: float = 10.0
    necrotic_thresh: float = 2.0
    tumor_growth_start_mcs: int = 1000
    vascular_vegf_activation_threshold: float = 0.2
    inactive_neighbor_area_limit: float = 200
    active_neighbor_area_limit: float = 150

    # Paper-derived mechanics
    tumor_target_volume: float = 25
    tumor_lambda_volume: float = 20.0
    necrotic_target_volume: float = 5.0
    tumor_target_surface: float = 2.0
    tumor_lambda_surface: float = 2.0
    vascular_target_volume: float = 55.0
    vascular_lambda_volume: float = 15.0
    vascular_target_surface: float = 120.0
    vascular_lambda_surface: float = 3.0
    necrotic_volume_loss_rate: float = 0.5

    # Vascular activation: VEGF2 threshold for switching Inactive -> Active
    vascular_activation_vegf2_threshold: float = 0.3
    # Vascular deactivation: VEGF2 threshold for switching Active -> Inactive
    vascular_deactivation_vegf2_threshold: float = 0.05

    # Paper-derived growth laws
    tumor_growth_volume_rate: float = 0.1
    tumor_growth_surface_rate: float = 0.1
    tumor_growth_denominator: float = 15.0
    vascular_growth_volume_rate: float = 0.2
    vascular_growth_surface_rate: float = 0.3
    vascular_growth_denominator: float = 1.0

    # Paper-derived mitosis thresholds
    tumor_doubling_volume: float = 50
    vascular_doubling_volume: float = 80.0

    # Field clamps for Python-side monitoring / safety
    field_max_value: float = 100.0
    report_frequency: int = 50

    # Python monitoring toggles
    monitor_frequency: int = 50
    monitor_to_console: bool = True
    monitor_to_csv: bool = True
    monitor_include_growth_rates: bool = True
    monitor_include_field_means: bool = True
    monitor_include_vascular_metrics: bool = True
    monitor_output_dir: str = "results"
    monitor_output_filename: str = "angiogenesis_metrics.csv"


BASE_CONFIG = ModelConfig()

PRESETS = {
    "paper_mvp": replace(
        BASE_CONFIG,
        preset_name="paper_mvp",
        enable_initial_sprouts=False,
        enable_type_switching=False,
        enable_tumor_growth=False,
        enable_vascular_growth=False,
        enable_mitosis=False,
        monitor_include_growth_rates=False,
        monitor_include_vascular_metrics=False,
    ),
    "paper_hypoxia": replace(
        BASE_CONFIG,
        preset_name="paper_hypoxia",
        enable_initial_sprouts=False,
        enable_type_switching=True,
        enable_tumor_growth=True,
        enable_vascular_growth=False,
        enable_mitosis=False,
        monitor_include_vascular_metrics=False,
    ),
    "paper_sprouting": replace(
        BASE_CONFIG,
        preset_name="paper_sprouting",
        enable_initial_sprouts=True,
        enable_type_switching=True,
        enable_tumor_growth=True,
        enable_vascular_growth=True,
        enable_mitosis=False,
    ),
    "tumor_oxygen_only": replace(
        BASE_CONFIG,
        preset_name="tumor_oxygen_only",
        enable_initial_sprouts=False,
        enable_type_switching=True,
        enable_tumor_growth=True,
        enable_vascular_growth=False,
        enable_mitosis=True,
        monitor_include_vascular_metrics=False,
        enable_timescale_separation=False,
    ),
    "paper_full": replace(
        BASE_CONFIG,
        preset_name="paper_full",
    ),
}

# Change this string to step through the model from simple to complex.
SELECTED_PRESET = "tumor_oxygen_only"
CONFIG = PRESETS[SELECTED_PRESET]
