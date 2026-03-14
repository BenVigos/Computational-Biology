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
    enable_hif1a_network: bool = True

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
    tip_cell_center_y_fractions: tuple[float, ...] = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

    # Bottom-vessel sprouts: InactiveNeovascular cells attached to the bottom boundary vessel
    enable_bottom_vessel_sprouts: bool = True
    bottom_sprout_x_fractions: tuple[float, ...] = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
    bottom_sprout_width_fraction: float = 0.04
    bottom_sprout_height_fraction: float = 0.04

    # Initial tumor geometry: a small circular cluster of cells
    tumor_center_x_fraction: float = 0.50
    tumor_center_y_fraction: float = 0.50
    tumor_radius_fraction: float = 0.05
    tumor_seed_size_fraction: float = 0.02

    # Paper-derived phenotype thresholds and timing
    area_thresh: float = 100.0
    nutrient_thresh: float = 20.0
    necrotic_thresh: float = 10.0
    tumor_growth_start_mcs: int = 500
    vascular_vegf_activation_threshold: float = 0.3
    inactive_neighbor_area_limit: float = 30
    active_neighbor_area_limit: float = 20

    # Paper-derived mechanics
    tumor_target_volume: float = 22.5
    tumor_lambda_volume: float = 8.0
    necrotic_target_volume: float = 0.0
    tumor_target_surface: float = 2.0
    tumor_lambda_surface: float = 2.0
    vascular_target_volume: float = 40.0
    vascular_lambda_volume: float = 15.0
    vascular_target_surface: float = 120.0
    vascular_lambda_surface: float = 3.0
    necrotic_volume_loss_rate: float = 0.5

    # Vascular activation: VEGF2 threshold for switching Inactive -> Active
    vascular_activation_vegf2_threshold: float = 0.3
    # Vascular deactivation: VEGF2 threshold for switching Active -> Inactive
    vascular_deactivation_vegf2_threshold: float = 0.05

    # Paper-derived growth laws
    tumor_growth_volume_rate: float = 0.15
    tumor_growth_surface_rate: float = 0.1
    tumor_growth_denominator: float = 10.0
    hypoxic_growth_volume_rate: float = tumor_growth_volume_rate**0.4
    hypoxic_growth_denominator: float = 10.0
    vascular_growth_volume_rate: float = 0.2
    vascular_growth_surface_rate: float = 0.3
    vascular_growth_denominator: float = 0.5

    # Oxygen -> HIF-1a -> VEGF intracellular signaling proxy
    hif1a_initial_value: float = 0.0
    hif1a_max_value: float = 1.0
    hif1a_hypoxia_halfmax_oxygen: float = 20.0
    hif1a_stabilization_rate: float = 0.04
    hif1a_degradation_rate: float = 0.05
    vegf_drive_basal: float = 0.02
    vegf_drive_max: float = 1.0
    vegf_drive_hill_k: float = 0.35
    vegf_drive_hill_n: float = 2.0
    hif1a_to_vegf2_weight: float = 0.6

    # Paper-derived mitosis thresholds
    tumor_doubling_volume: float = 50
    vascular_doubling_volume: float = 80.0

    # Field clamps for Python-side monitoring / safety
    field_max_value: float = 100.0
    report_frequency: int = 50

    # Python monitoring toggles
    monitor_frequency: int = 10
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
    "vegf_tuning": replace(
        BASE_CONFIG,
        preset_name="tumor_oxygen_only",
        enable_initial_sprouts=False,
        enable_bottom_vessel_sprouts=False,
        enable_type_switching=True,
        enable_tumor_growth=True,
        enable_vascular_growth=False,
        enable_mitosis=True,
        monitor_include_vascular_metrics=False,
        enable_timescale_separation=False,
        # nutrient_thresh=100,
        # necrotic_thresh=0.0,
        tumor_growth_start_mcs=300,
        tumor_radius_fraction=0.16,
    ),

    "branching_tuning": replace(
        BASE_CONFIG,
        preset_name="branching_tuning",
        enable_initial_sprouts=False,
        enable_bottom_vessel_sprouts=True,
        enable_type_switching=True,
        enable_tumor_growth=True,
        enable_vascular_growth=True,
        enable_mitosis=True,
        monitor_include_vascular_metrics=True,
        enable_timescale_separation=False,
        vascular_activation_vegf2_threshold=0.7,
        vascular_deactivation_vegf2_threshold=0.4,
        nutrient_thresh=22.5,
        necrotic_thresh=10,
        tumor_growth_start_mcs=300,
        tumor_radius_fraction=0.16,
        vascular_growth_volume_rate=0.5,
    ),

    "paper_full": replace(
        BASE_CONFIG,
        preset_name="paper_full",
    ),
}

# Change this string to step through the model from simple to complex.
SELECTED_PRESET = "branching_tuning"
CONFIG = PRESETS[SELECTED_PRESET]
