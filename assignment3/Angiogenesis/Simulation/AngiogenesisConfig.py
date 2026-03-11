from dataclasses import dataclass


@dataclass(frozen=True)
class ModelConfig:
    preset_name: str = "paper_2d"

    # Runtime toggles
    enable_type_switching: bool = True
    enable_tumor_growth: bool = True
    enable_vascular_growth: bool = True
    enable_mitosis: bool = True
    enable_python_monitoring: bool = True

    # 2D lattice / initialization
    lattice_x: int = 180
    lattice_y: int = 180
    steps: int = 5000
    vessel_x_min: int = 18
    vessel_x_max: int = 28
    vessel_y_min: int = 24
    vessel_y_max: int = 156
    tip_cell_x_min: int = 28
    tip_cell_x_max: int = 36
    tip_cell_y_values: tuple[int, ...] = (42, 90, 138)
    tip_cell_height: int = 10
    tumor_x_min: int = 72
    tumor_x_max: int = 108
    tumor_y_min: int = 72
    tumor_y_max: int = 108
    tumor_seed_size: int = 6

    # Paper-derived phenotype thresholds and timing
    area_thresh: float = 1.0
    nutrient_thresh: float = 5.0
    necrotic_thresh: float = 1.0
    tumor_growth_start_mcs: int = 100
    vascular_vegf_activation_threshold: float = 0.5
    inactive_neighbor_area_limit: float = 70.0
    active_neighbor_area_limit: float = 50.0

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
    tumor_growth_volume_rate: float = 0.04
    tumor_growth_surface_rate: float = 0.12
    tumor_growth_denominator: float = 10.0
    vascular_growth_volume_rate: float = 0.06
    vascular_growth_surface_rate: float = 0.15
    vascular_growth_denominator: float = 0.5

    # Paper-derived mitosis thresholds
    tumor_doubling_volume: float = 54.0
    vascular_doubling_volume: float = 80.0

    # Field clamps for Python-side monitoring / safety
    field_max_value: float = 100.0
    report_frequency: int = 100

    # Python monitoring toggles
    monitor_frequency: int = 50
    monitor_to_console: bool = True
    monitor_to_csv: bool = True
    monitor_include_growth_rates: bool = True
    monitor_include_field_means: bool = True
    monitor_include_vascular_metrics: bool = True
    monitor_output_dir: str = "results"
    monitor_output_filename: str = "angiogenesis_metrics.csv"


CONFIG = ModelConfig()
