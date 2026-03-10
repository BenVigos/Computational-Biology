from dataclasses import dataclass, replace


@dataclass(frozen=True)
class ModelConfig:
    preset_name: str = "mvp"

    # Feature toggles
    enable_oxygen: bool = True
    enable_oxygen_uptake: bool = True
    enable_hypoxia_switch: bool = True
    enable_hif: bool = False
    enable_vegf_secretion: bool = True
    enable_tumor_growth: bool = True
    enable_endothelial_growth: bool = False
    enable_chemotaxis: bool = False
    enable_mitosis: bool = False
    enable_apoptosis: bool = False

    # Lattice / initialization
    lattice_x: int = 160
    lattice_y: int = 160
    tumor_x_min: int = 70
    tumor_x_max: int = 98
    tumor_y_min: int = 66
    tumor_y_max: int = 94
    tumor_seed_size: int = 7
    vessel_x_min: int = 50
    vessel_x_max: int = 60
    vessel_y_min: int = 20
    vessel_y_max: int = 140
    tip_cell_x_min: int = 60
    tip_cell_x_max: int = 70
    tip_cell_y_values: tuple[int, ...] = (52, 76, 100)
    tip_cell_height: int = 8

    # Mechanics
    default_target_volume: float = 25.0
    default_lambda_volume: float = 20
    blood_vessel_lambda_volume: float = 8.0
    endothelial_lambda_volume: float = 3.0

    # Oxygen dynamics
    oxygen_supply_rate: float = 0.08
    oxygen_uptake_rate: float = 0.015
    oxygen_hypoxia_threshold: float = 0.12
    oxygen_death_threshold: float = 0.2
    oxygen_mm_constant: float = 0.10

    # Tumor growth
    tumor_growth_gm: float = 1.0
    apoptosis_rate: float = 0.1

    # HIF-1 alpha dynamics
    hif_alpha: float = 0.02
    hif_beta: float = 0.40
    hif_to_vegf_gain: float = 1.0

    # VEGF / endothelial response
    vegf_secretion_rate: float = 0.03
    vegf_activation_threshold: float = 0.10
    vegf_half_max_scale: float = 2.0
    endothelial_growth_gv: float = 0.10
    chemotaxis_lambda: float = 120.0

    # Division thresholds
    tumor_division_volume: float = 50
    endothelial_division_volume: float = 70.0

    # Diagnostics / clamps
    field_max_value: float = 5.0
    report_frequency: int = 100

    # Python-side monitoring toggles
    # These control optional time-series metrics computed in Python, such as
    # average tumor volume growth rate, hypoxic fraction, and local field means.
    enable_python_monitoring: bool = True
    monitor_frequency: int = 50
    monitor_to_console: bool = True
    monitor_to_csv: bool = True
    monitor_include_growth_rates: bool = True
    monitor_include_field_means: bool = True
    monitor_include_endothelial_metrics: bool = True
    monitor_output_dir: str = "results"
    monitor_output_filename: str = "angiogenesis_metrics.csv"


BASE_CONFIG = ModelConfig()

PRESETS = {
    "mvp": BASE_CONFIG,
    "hypoxia_hif": replace(BASE_CONFIG, preset_name="hypoxia_hif", enable_hif=True),
    "sprouting": replace(
        BASE_CONFIG,
        preset_name="sprouting",
        enable_hif=True,
        enable_endothelial_growth=True,
        enable_chemotaxis=True,
    ),
    "full": replace(
        BASE_CONFIG,
        preset_name="full",
        enable_hif=True,
        enable_endothelial_growth=True,
        enable_chemotaxis=True,
        enable_mitosis=True,
        enable_apoptosis=True,
    ),
}

# Change this one string to grow the model in stages.
SELECTED_PRESET = "full"
CONFIG = PRESETS[SELECTED_PRESET]
