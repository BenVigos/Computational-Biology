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

    # 2D lattice / initialization
    lattice_x: int = 100
    lattice_y: int = 100
    steps: int = 3000
    vessel_x_min: int = 0
    vessel_x_max: int = 8
    vessel_y_min: int = 0
    vessel_y_max: int = 100
    tip_cell_x_min: int = vessel_x_max + 1
    tip_cell_x_max: int = tip_cell_x_min + 5 
    tip_cell_y_values: tuple[int, ...] = (16, 22, 28, 34, 40, 48, 54, 60, 66, 74, 80, 86)
    tip_cell_height: int = 6
    tumor_x_min: int = 45
    tumor_x_max: int = 55
    tumor_y_min: int = 45
    tumor_y_max: int = 55
    tumor_seed_size: int = 1

    # Paper-derived phenotype thresholds and timing
    area_thresh: float = 1.0
    nutrient_thresh: float = 5.0
    necrotic_thresh: float = 1.0
    tumor_growth_start_mcs: int = 5
    vascular_vegf_activation_threshold: float = 0.5
    inactive_neighbor_area_limit: float = 7000.0
    active_neighbor_area_limit: float = 5000.0

    # Paper-derived mechanics
    tumor_target_volume: float = 33.0
    tumor_lambda_volume: float = 10.0
    tumor_target_surface: float = 90.0
    tumor_lambda_surface: float = 2.0
    vascular_target_volume: float = 60.0
    vascular_lambda_volume: float = 13.0
    vascular_target_surface: float = 150.0
    vascular_lambda_surface: float = 3.0
    necrotic_volume_loss_rate: float = 0.5

    # Paper-derived growth laws
    tumor_growth_volume_rate: float = 0.02
    tumor_growth_surface_rate: float = 0.06
    tumor_growth_denominator: float = 10.0
    vascular_growth_volume_rate: float = 0.06
    vascular_growth_surface_rate: float = 0.15
    vascular_growth_denominator: float = 0.5

    # Paper-derived mitosis thresholds
    tumor_doubling_volume: float = 80.0
    vascular_doubling_volume: float = 120.0

    # Field clamps for Python-side monitoring / safety
    field_max_value: float = 100.0
    report_frequency: int = 100

    # Python monitoring toggles
    monitor_frequency: int = 1
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
    "paper_full": replace(
        BASE_CONFIG,
        preset_name="paper_full",
    ),
}

# Change this string to step through the model from simple to complex.
SELECTED_PRESET = "paper_sprouting"
CONFIG = PRESETS[SELECTED_PRESET]
