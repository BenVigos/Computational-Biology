# Paper-aligned 2D angiogenesis steppables.
# The mechanics mirror `Paper files/Simulation/TumorVasc3DSteppables.py`
# while keeping the lattice 2D and the initial geometry programmatic.

import csv
import math
import pathlib

from cc3d.core.PySteppables import *

from AngiogenesisConfig import CONFIG


def _timescale_cycle_length():
    diffusion_mcs = max(0, int(CONFIG.diffusion_relaxation_mcs))
    growth_mcs = max(0, int(CONFIG.growth_window_mcs))
    total = diffusion_mcs + growth_mcs
    return total if total > 0 else 1


def _is_growth_window(mcs):
    if not CONFIG.enable_timescale_separation:
        return True

    diffusion_mcs = max(0, int(CONFIG.diffusion_relaxation_mcs))
    growth_mcs = max(0, int(CONFIG.growth_window_mcs))
    if growth_mcs <= 0:
        return False

    phase_index = (int(mcs) + int(CONFIG.timescale_cycle_offset_mcs)) % _timescale_cycle_length()
    return phase_index >= diffusion_mcs


def _timescale_phase_label(mcs):
    return "growth" if _is_growth_window(mcs) else "diffusion"


class BaseModelSteppable(SteppableBasePy):
    NORMAL: int
    HYPOXIC: int
    NECROTIC: int
    ACTIVENEOVASCULAR: int
    VASCULAR: int
    INACTIVENEOVASCULAR: int

    def __init__(self, frequency=1):
        super().__init__(frequency)

    def _clamp(self, value, minimum=0.0, maximum=None):
        if maximum is not None and value > maximum:
            return maximum
        if value < minimum:
            return minimum
        return value

    def _cell_com_indices(self, cell):
        volume = max(float(getattr(cell, "volume", 0.0)), 1e-6)
        if hasattr(cell, "xCOM"):
            x = int(round(float(cell.xCOM)))
            y = int(round(float(cell.yCOM)))
            z = int(round(float(cell.zCOM)))
        else:
            x = int(round(float(cell.xCM) / volume))
            y = int(round(float(cell.yCM) / volume))
            z = int(round(float(cell.zCM) / volume))

        dim_x = int(getattr(self.dim, "x"))
        dim_y = int(getattr(self.dim, "y"))
        dim_z = int(getattr(self.dim, "z"))
        return (
            min(max(x, 0), dim_x - 1),
            min(max(y, 0), dim_y - 1),
            min(max(z, 0), dim_z - 1),
        )

    def _field_at_com(self, field, cell):
        x, y, z = self._cell_com_indices(cell)
        return float(field[x, y, z])

    def _mean(self, values):
        return sum(values) / len(values) if values else 0.0

    def _cell_com_position(self, cell):
        x, y, _ = self._cell_com_indices(cell)
        return float(x), float(y)

    def _weighted_centroid(self, cells):
        if not cells:
            return None

        total_weight = 0.0
        sum_x = 0.0
        sum_y = 0.0
        for cell in cells:
            weight = max(float(getattr(cell, "volume", 0.0)), 1.0)
            x, y = self._cell_com_position(cell)
            sum_x += weight * x
            sum_y += weight * y
            total_weight += weight

        if total_weight <= 0.0:
            return None
        return sum_x / total_weight, sum_y / total_weight

    def _euclidean_distance(self, point_a, point_b):
        if point_a is None or point_b is None:
            return 0.0
        dx = float(point_a[0]) - float(point_b[0])
        dy = float(point_a[1]) - float(point_b[1])
        return math.sqrt(dx * dx + dy * dy)

    def _min_distance_to_points(self, source_point, target_points):
        if source_point is None or not target_points:
            return 0.0
        return min(self._euclidean_distance(source_point, point) for point in target_points)

    def _radial_distances(self, cells, origin):
        if origin is None:
            return []
        return [self._euclidean_distance(self._cell_com_position(cell), origin) for cell in cells]

    def _vascular_type_ids(self):
        return (self.ACTIVENEOVASCULAR, self.VASCULAR, self.INACTIVENEOVASCULAR)

    def _tumor_type_ids(self):
        return (self.NORMAL, self.HYPOXIC, self.NECROTIC)

    def _growing_tumor_type_ids(self):
        return (self.NORMAL, self.HYPOXIC)

    def _dividing_type_ids(self):
        return (self.NORMAL, self.HYPOXIC, self.ACTIVENEOVASCULAR, self.INACTIVENEOVASCULAR)

    def _apply_paper_mechanics(self, cell):
        if cell.type in self._vascular_type_ids():
            cell.targetVolume = CONFIG.vascular_target_volume
            cell.lambdaVolume = CONFIG.vascular_lambda_volume

        elif cell.type == self.NECROTIC:
            cell.targetVolume = CONFIG.necrotic_target_volume
            cell.lambdaVolume = 8.0
        else:
            cell.targetVolume = CONFIG.tumor_target_volume
            cell.lambdaVolume = CONFIG.tumor_lambda_volume

    def _fraction_to_index(self, fraction, size, *, inclusive_upper=False):
        fraction = self._clamp(float(fraction), minimum=0.0, maximum=1.0)
        upper = size - 1 if inclusive_upper else size
        return int(round(fraction * upper))

    def _fraction_to_span(self, start_fraction, end_fraction, size, *, min_size=1):
        start = self._fraction_to_index(start_fraction, size)
        end = self._fraction_to_index(end_fraction, size)
        start = min(max(start, 0), size)
        end = min(max(end, 0), size)
        if end < start:
            start, end = end, start
        if end == start:
            end = min(size, start + max(min_size, 1))
            start = max(0, end - max(min_size, 1))
        return start, end

    def _fraction_to_length(self, length_fraction, size, *, min_size=1):
        length = int(round(self._clamp(float(length_fraction), minimum=0.0, maximum=1.0) * size))
        return min(size, max(min_size, length))

    def _is_growth_window(self, mcs):
        return _is_growth_window(mcs)

    def _timescale_phase_label(self, mcs):
        return _timescale_phase_label(mcs)


class ConstraintInitializerSteppable(BaseModelSteppable):
    def start(self):
        dim_x = int(getattr(self.dim, "x"))
        dim_y = int(getattr(self.dim, "y"))
        seeded_tumor_cells = 0

        vessel_thickness = self._fraction_to_length(
            CONFIG.vessel_boundary_thickness_fraction,
            min(dim_x, dim_y),
            min_size=1,
        )
        vessel_margin_x = self._fraction_to_index(CONFIG.vessel_boundary_margin_fraction, dim_x)
        vessel_margin_y = self._fraction_to_index(CONFIG.vessel_boundary_margin_fraction, dim_y)

        if CONFIG.enable_initial_vascular_strip:
            y0, y1 = vessel_margin_y, max(vessel_margin_y + 1, dim_y - vessel_margin_y)
            x0, x1 = vessel_margin_x, max(vessel_margin_x + 1, dim_x - vessel_margin_x)

            left_vessel = self.new_cell(self.VASCULAR)
            self.cell_field[0:vessel_thickness, y0:y1, 0] = left_vessel
            self._apply_paper_mechanics(left_vessel)

            right_vessel = self.new_cell(self.VASCULAR)
            self.cell_field[max(0, dim_x - vessel_thickness):dim_x, y0:y1, 0] = right_vessel
            self._apply_paper_mechanics(right_vessel)

            top_vessel = self.new_cell(self.VASCULAR)
            self.cell_field[x0:x1, 0:vessel_thickness, 0] = top_vessel
            self._apply_paper_mechanics(top_vessel)

            bottom_vessel = self.new_cell(self.VASCULAR)
            self.cell_field[x0:x1, max(0, dim_y - vessel_thickness):dim_y, 0] = bottom_vessel
            self._apply_paper_mechanics(bottom_vessel)

            if CONFIG.enable_top_vessel_sprouts:
                sprout_w = self._fraction_to_length(
                    CONFIG.bottom_sprout_width_fraction, dim_x, min_size=1
                )
                sprout_h = self._fraction_to_length(
                    CONFIG.bottom_sprout_height_fraction, dim_y, min_size=1
                )
                # sprouts sit immediately above (inside) the bottom vessel band
                sprout_y_max = max(0, dim_y - vessel_thickness)
                sprout_y_min = max(0, sprout_y_max - sprout_h)

                for x_center_fraction in CONFIG.bottom_sprout_x_fractions:
                    sx_center = self._fraction_to_index(
                        x_center_fraction, dim_x, inclusive_upper=True
                    )
                    sx_min = max(0, sx_center - sprout_w // 2)
                    sx_max = min(dim_x, sx_min + sprout_w)
                    sx_min = max(0, sx_max - sprout_w)  # re-clamp after upper bound

                    sprout = self.new_cell(self.INACTIVENEOVASCULAR)
                    self.cell_field[sx_min:sx_max, sprout_y_min:sprout_y_max, 0] = sprout
                    self._apply_paper_mechanics(sprout)

        if CONFIG.enable_initial_sprouts:
            tip_width = self._fraction_to_length(CONFIG.tip_cell_width_fraction, dim_x, min_size=1)
            tip_height = self._fraction_to_length(CONFIG.tip_cell_height_fraction, dim_y, min_size=1)
            tip_x_min = min(
                max(0, vessel_thickness + self._fraction_to_index(CONFIG.tip_cell_x_offset_fraction, dim_x)),
                max(0, dim_x - 1),
            )
            tip_x_max = min(dim_x, tip_x_min + tip_width)

            for y_center_fraction in CONFIG.tip_cell_center_y_fractions:
                y_center = self._fraction_to_index(y_center_fraction, dim_y, inclusive_upper=True)
                y_min = max(0, y_center - tip_height // 2)
                y_max = min(dim_y, y_min + tip_height)
                y_min = max(0, y_max - tip_height)
                tip = self.new_cell(self.INACTIVENEOVASCULAR)
                self.cell_field[tip_x_min:tip_x_max, y_min:y_max, 0] = tip
                self._apply_paper_mechanics(tip)

        if CONFIG.enable_initial_tumor_mass:
            tumor_center_x = self._fraction_to_index(
                CONFIG.tumor_center_x_fraction,
                dim_x,
                inclusive_upper=True,
            )
            tumor_center_y = self._fraction_to_index(
                CONFIG.tumor_center_y_fraction,
                dim_y,
                inclusive_upper=True,
            )
            tumor_radius = self._fraction_to_length(
                CONFIG.tumor_radius_fraction,
                min(dim_x, dim_y),
                min_size=1,
            )
            tumor_seed_size = self._fraction_to_length(
                CONFIG.tumor_seed_size_fraction,
                min(dim_x, dim_y),
                min_size=1,
            )
            radius_squared = float(tumor_radius * tumor_radius)

            x_start = max(0, tumor_center_x - tumor_radius)
            x_stop = min(dim_x, tumor_center_x + tumor_radius + 1)
            y_start = max(0, tumor_center_y - tumor_radius)
            y_stop = min(dim_y, tumor_center_y + tumor_radius + 1)

            for x_min in range(x_start, x_stop, tumor_seed_size):
                for y_min in range(y_start, y_stop, tumor_seed_size):
                    x_max = min(x_min + tumor_seed_size, dim_x)
                    y_max = min(y_min + tumor_seed_size, dim_y)
                    patch_center_x = 0.5 * (x_min + x_max - 1)
                    patch_center_y = 0.5 * (y_min + y_max - 1)

                    if (
                        (patch_center_x - tumor_center_x) ** 2
                        + (patch_center_y - tumor_center_y) ** 2
                    ) > radius_squared:
                        continue

                    tumor = self.new_cell(self.NORMAL)
                    self.cell_field[x_min:x_max, y_min:y_max, 0] = tumor
                    self._apply_paper_mechanics(tumor)
                    seeded_tumor_cells += 1

        print(
            f"[Angiogenesis] Initialised preset '{CONFIG.preset_name}' on a {dim_x}x{dim_y} grid "
            f"with {seeded_tumor_cells} initial tumor cells."
        )


class GrowthSteppable(BaseModelSteppable):
    def start(self):
        for cell in self.cellList:
            self._apply_paper_mechanics(cell)
        if CONFIG.enable_hif1a_network:
            for cell in self.cell_list_by_type(*self._tumor_type_ids()):
                cell.dict["hif1a"] = float(CONFIG.hif1a_initial_value)
                cell.dict["vegf_drive"] = float(CONFIG.vegf_drive_basal)

    def _oxygen_to_hypoxia_signal(self, local_oxygen):
        # Low oxygen yields a larger hypoxia signal in [0, 1].
        oxygen = max(0.0, float(local_oxygen))
        k = max(1e-6, float(CONFIG.hif1a_hypoxia_halfmax_oxygen))
        return k / (k + oxygen)

    def _hill_activation(self, value, k, n):
        value = max(0.0, float(value))
        k = max(1e-6, float(k))
        n = max(1.0, float(n))
        value_power = value ** n
        return value_power / (k ** n + value_power)

    def _update_intracellular_network(self, cell, local_oxygen):
        if not CONFIG.enable_hif1a_network or cell.type not in self._tumor_type_ids():
            return

        current_hif1a = float(cell.dict.get("hif1a", CONFIG.hif1a_initial_value))
        hypoxia_signal = self._oxygen_to_hypoxia_signal(local_oxygen)

        stabilization = float(CONFIG.hif1a_stabilization_rate) * hypoxia_signal
        degradation = float(CONFIG.hif1a_degradation_rate) * current_hif1a
        updated_hif1a = self._clamp(
            current_hif1a + stabilization - degradation,
            minimum=0.0,
            maximum=float(CONFIG.hif1a_max_value),
        )

        vegf_hif_fraction = self._hill_activation(
            updated_hif1a,
            CONFIG.vegf_drive_hill_k,
            CONFIG.vegf_drive_hill_n,
        )
        updated_vegf_drive = self._clamp(
            float(CONFIG.vegf_drive_basal)
            + (float(CONFIG.vegf_drive_max) - float(CONFIG.vegf_drive_basal)) * vegf_hif_fraction,
            minimum=0.0,
            maximum=float(CONFIG.vegf_drive_max),
        )

        cell.dict["hif1a"] = updated_hif1a
        cell.dict["vegf_drive"] = updated_vegf_drive

    def _mean_tumor_vegf_drive(self):
        tumor_cells = list(self.cell_list_by_type(*self._tumor_type_ids()))
        if not tumor_cells:
            return 0.0
        return self._mean([float(cell.dict.get("vegf_drive", 0.0)) for cell in tumor_cells])

    def _effective_vegf2_signal(self, local_vegf2):
        if not CONFIG.enable_hif1a_network:
            return float(local_vegf2)
        return max(
            0.0,
            float(local_vegf2) + float(CONFIG.hif1a_to_vegf2_weight) * self._mean_tumor_vegf_drive(),
        )

    def _neovascular_neighbor_area(self, cell):
        """Shared surface with other *neovascular* cells only (Active + Inactive).
        The frozen Vascular strip is excluded so that touching the parent vessel
        does not block growth."""
        total_area = 0.0
        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
            if neighbor and neighbor.type in (self.ACTIVENEOVASCULAR, self.INACTIVENEOVASCULAR):
                total_area += float(common_surface_area)
        return total_area

    def _vascular_growth_step(self, cell, local_vegf2, neighbor_area_limit):
        if local_vegf2 <= CONFIG.vascular_vegf_activation_threshold:
            return

        if self._neovascular_neighbor_area(cell) >= neighbor_area_limit:
            return

        growth_fraction = local_vegf2 / (CONFIG.vascular_growth_denominator + local_vegf2)
        cell.targetVolume += CONFIG.vascular_growth_volume_rate * growth_fraction

    def _vascular_activation_step(self, cell, local_vegf2):
        """Switch InactiveNeovascular -> ActiveNeovascular when VEGF2 is high,
        and ActiveNeovascular -> InactiveNeovascular when VEGF2 drops."""
        if not CONFIG.enable_type_switching:
            return

        if cell.type == self.INACTIVENEOVASCULAR:
            if local_vegf2 > CONFIG.vascular_activation_vegf2_threshold:
                cell.type = self.ACTIVENEOVASCULAR
        elif cell.type == self.ACTIVENEOVASCULAR:
            if local_vegf2 < CONFIG.vascular_deactivation_vegf2_threshold:
                cell.type = self.INACTIVENEOVASCULAR

    def _tumor_growth_step(self, cell, local_oxygen, mcs):
        if CONFIG.enable_type_switching:
            if local_oxygen < CONFIG.nutrient_thresh and mcs > CONFIG.tumor_growth_start_mcs:
                cell.type = self.HYPOXIC

            if local_oxygen < CONFIG.necrotic_thresh and mcs > CONFIG.tumor_growth_start_mcs:
                cell.type = self.NECROTIC

        if CONFIG.enable_tumor_growth and mcs > CONFIG.tumor_growth_start_mcs and cell.type == self.NORMAL:
            cell.targetVolume += (
                CONFIG.tumor_growth_volume_rate
                * local_oxygen
                / (CONFIG.tumor_growth_denominator + local_oxygen)
            )

        if CONFIG.enable_tumor_growth and mcs > CONFIG.tumor_growth_start_mcs and cell.type == self.HYPOXIC:
            cell.targetVolume += (
                CONFIG.hypoxic_growth_volume_rate
                * local_oxygen
                / (CONFIG.hypoxic_growth_denominator + local_oxygen)
            )
        if CONFIG.enable_tumor_growth and mcs > CONFIG.tumor_growth_start_mcs and cell.type== self.NECROTIC:
            cell.targetVolume = 0.0
            cell.lambdaVolume = 8.0

    def _hypoxic_adjustments(self, cell, local_oxygen, mcs):
        if not CONFIG.enable_type_switching:
            return

        if local_oxygen < CONFIG.necrotic_thresh and mcs > CONFIG.tumor_growth_start_mcs:
            cell.type = self.NECROTIC
        elif local_oxygen > CONFIG.nutrient_thresh:
            cell.type = self.NORMAL

    def _necrotic_adjustments(self, cell):
        # cell.targetVolume = self._clamp(
        #     float(cell.targetVolume) - CONFIG.necrotic_volume_loss_rate,
        #     minimum=0.0,
        # )
        cell.targetVolume = 0.0
        cell.lambdaVolume = 8.0
        if cell.targetVolume < 5.0:
            cell.targetVolume = 0.0

    def step(self, mcs):
        if not self._is_growth_window(mcs):
            return

        field_neo_vasc = self.field.VEGF2
        field_malig = self.field.Oxygen

        for cell in self.cellList:
            if cell.type in self._tumor_type_ids():
                self._update_intracellular_network(cell, self._field_at_com(field_malig, cell))

            if cell.type == self.INACTIVENEOVASCULAR and CONFIG.enable_vascular_growth:
                local_vegf2 = self._field_at_com(field_neo_vasc, cell)
                effective_vegf2 = self._effective_vegf2_signal(local_vegf2)
                self._vascular_activation_step(cell, effective_vegf2)
                self._vascular_growth_step(
                    cell,
                    effective_vegf2,
                    CONFIG.inactive_neighbor_area_limit,
                )

            elif cell.type == self.ACTIVENEOVASCULAR and CONFIG.enable_vascular_growth:
                local_vegf2 = self._field_at_com(field_neo_vasc, cell)
                effective_vegf2 = self._effective_vegf2_signal(local_vegf2)
                self._vascular_activation_step(cell, effective_vegf2)
                self._vascular_growth_step(
                    cell,
                    effective_vegf2,
                    CONFIG.active_neighbor_area_limit,
                )

            elif cell.type in self._growing_tumor_type_ids():
                self._tumor_growth_step(cell, self._field_at_com(field_malig, cell), mcs)

            if cell.type == self.HYPOXIC:
                self._hypoxic_adjustments(cell, self._field_at_com(field_malig, cell), mcs)

            if cell.type == self.NECROTIC:
                self._necrotic_adjustments(cell)


class MitosisSteppable(MitosisSteppableBase):
    NORMAL: int
    HYPOXIC: int
    NECROTIC: int
    ACTIVENEOVASCULAR: int
    VASCULAR: int
    INACTIVENEOVASCULAR: int

    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.doubling_volume_dict = {
            1: CONFIG.tumor_doubling_volume,
            2: CONFIG.tumor_doubling_volume,
            4: CONFIG.vascular_doubling_volume,
            6: CONFIG.vascular_doubling_volume,
        }

    def step(self, mcs):
        if not CONFIG.enable_mitosis:
            return
        if not _is_growth_window(mcs):
            return

        cells_to_divide = []
        for cell in self.cell_list_by_type(*self.doubling_volume_dict.keys()):
            doubling_volume = self.doubling_volume_dict.get(cell.type)
            if doubling_volume is not None and cell.volume > doubling_volume:
                cells_to_divide.append(cell)

        for cell in cells_to_divide:
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        parent_cell = getattr(self, "parent_cell", None) or getattr(self, "parentCell")
        child_cell = getattr(self, "child_cell", None) or getattr(self, "childCell")

        if parent_cell.type in (self.NORMAL, self.HYPOXIC):
            child_cell.type = parent_cell.type
            child_cell.targetVolume = CONFIG.tumor_target_volume
            child_cell.lambdaVolume = CONFIG.tumor_lambda_volume
            parent_cell.targetVolume = CONFIG.tumor_target_volume
            parent_cell.lambdaVolume = CONFIG.tumor_lambda_volume
            if CONFIG.enable_hif1a_network:
                parent_hif1a = float(parent_cell.dict.get("hif1a", CONFIG.hif1a_initial_value))
                parent_vegf = float(parent_cell.dict.get("vegf_drive", CONFIG.vegf_drive_basal))
                parent_cell.dict["hif1a"] = parent_hif1a
                parent_cell.dict["vegf_drive"] = parent_vegf
                child_cell.dict["hif1a"] = parent_hif1a
                child_cell.dict["vegf_drive"] = parent_vegf


        if parent_cell.type in (self.INACTIVENEOVASCULAR, self.ACTIVENEOVASCULAR):
            child_cell.type = parent_cell.type
            child_cell.targetVolume = CONFIG.vascular_target_volume
            child_cell.lambdaVolume = CONFIG.vascular_lambda_volume
            parent_cell.targetVolume = CONFIG.vascular_target_volume
            parent_cell.lambdaVolume = CONFIG.vascular_lambda_volume


class MonitoringSteppable(BaseModelSteppable):
    NORMAL: int
    HYPOXIC: int
    NECROTIC: int
    ACTIVENEOVASCULAR: int
    VASCULAR: int
    INACTIVENEOVASCULAR: int

    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.metric_history = []
        self._previous_sample_mcs = None
        self._previous_tumor_volumes = {}
        self._previous_tumor_targets = {}
        self._previous_endothelial_volumes = {}
        self._previous_vascular_volumes = {}
        self._csv_path = None

    def start(self):
        if not CONFIG.enable_python_monitoring or not CONFIG.monitor_to_csv:
            return

        output_dir = pathlib.Path(__file__).resolve().parents[1] / CONFIG.monitor_output_dir
        output_dir.mkdir(parents=True, exist_ok=True)
        self._csv_path = output_dir / CONFIG.monitor_output_filename
        self._csv_path.write_text("", encoding="utf-8")
        print(f"[Angiogenesis] Monitoring CSV -> {self._csv_path}")

    def _cell_value_map(self, cells, *, attribute_name):
        value_map = {}
        for cell in cells:
            cell_id = getattr(cell, "id", None)
            if cell_id is not None:
                value_map[cell_id] = float(getattr(cell, attribute_name))
        return value_map

    def _mean_delta_rate(self, current_values, previous_values, delta_mcs):
        if delta_mcs <= 0 or not previous_values:
            return 0.0

        if not current_values:
            return 0.0

        # Population-average rate that includes division/new IDs and disappearing IDs
        # via the total-volume delta, then normalizes by current population size.
        delta_total = sum(current_values.values()) - sum(previous_values.values())
        return delta_total / (delta_mcs * max(len(current_values), 1))

    def _total_delta_rate(self, current_values, previous_values, delta_mcs):
        if delta_mcs <= 0 or not previous_values:
            return 0.0
        return (sum(current_values.values()) - sum(previous_values.values())) / delta_mcs

    def _write_metrics_to_csv(self, metrics):
        if self._csv_path is None:
            return

        write_header = not self._csv_path.exists() or self._csv_path.stat().st_size == 0
        with self._csv_path.open("a", newline="", encoding="utf-8") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=list(metrics.keys()))
            if write_header:
                writer.writeheader()
            writer.writerow(metrics)

    def _collect_metrics(self, mcs):
        normal_cells = list(self.cell_list_by_type(self.NORMAL))
        hypoxic_cells = list(self.cell_list_by_type(self.HYPOXIC))
        necrotic_cells = list(self.cell_list_by_type(self.NECROTIC))
        tumor_like_cells = normal_cells + hypoxic_cells + necrotic_cells
        active_neovascular_cells = list(self.cell_list_by_type(self.ACTIVENEOVASCULAR))
        vascular_cells = list(self.cell_list_by_type(self.VASCULAR))
        inactive_neovascular_cells = list(self.cell_list_by_type(self.INACTIVENEOVASCULAR))
        vascular_like_cells = active_neovascular_cells + vascular_cells + inactive_neovascular_cells
        parent_vessel_cells = list(self.cell_list_by_type(self.VASCULAR))
        neovascular_cells   = active_neovascular_cells + inactive_neovascular_cells
        endothelial_cells = active_neovascular_cells + inactive_neovascular_cells #duplicate
        all_vascular_cells  = neovascular_cells + parent_vessel_cells 
        
        delta_mcs = 0 if self._previous_sample_mcs is None else mcs - self._previous_sample_mcs
        normal_volumes = [float(cell.volume) for cell in normal_cells]
        hypoxic_volumes = [float(cell.volume) for cell in hypoxic_cells]
        necrotic_volumes = [float(cell.volume) for cell in necrotic_cells]
        active_neovascular_volumes = [float(cell.volume) for cell in active_neovascular_cells]
        inactive_neovascular_volumes = [float(cell.volume) for cell in inactive_neovascular_cells]
        tumor_volumes = normal_volumes + hypoxic_volumes + necrotic_volumes
        tumor_target_volumes = [float(cell.targetVolume) for cell in tumor_like_cells]
        neovascular_volumes = [float(cell.volume) for cell in neovascular_cells]
        endothelial_volumes = [float(cell.volume) for cell in endothelial_cells]
        parent_vessel_volumes = [float(cell.volume) for cell in vascular_cells]

        current_tumor_volumes = self._cell_value_map(tumor_like_cells, attribute_name="volume")
        current_tumor_targets = self._cell_value_map(tumor_like_cells, attribute_name="targetVolume")
        current_endothelial_volumes = self._cell_value_map(endothelial_cells, attribute_name="volume")
        current_vascular_volumes = self._cell_value_map(vascular_cells, attribute_name="volume")

        tumor_centroid = self._weighted_centroid(tumor_like_cells)
        hypoxic_core_cells = hypoxic_cells + necrotic_cells
        hypoxic_core_centroid = self._weighted_centroid(hypoxic_core_cells)
        tumor_positions = [self._cell_com_position(cell) for cell in tumor_like_cells]
        neovascular_positions = [self._cell_com_position(cell) for cell in neovascular_cells]
        tumor_radial_distances = self._radial_distances(tumor_like_cells, tumor_centroid)
        hypoxic_radial_distances = self._radial_distances(hypoxic_cells, tumor_centroid)
        necrotic_radial_distances = self._radial_distances(necrotic_cells, tumor_centroid)
        tumor_effective_radius = math.sqrt(sum(tumor_volumes) / math.pi) if tumor_volumes else 0.0
        tumor_min_x = min((point[0] for point in tumor_positions), default=0.0)
        tumor_max_x = max((point[0] for point in tumor_positions), default=0.0)
        tumor_min_y = min((point[1] for point in tumor_positions), default=0.0)
        tumor_max_y = max((point[1] for point in tumor_positions), default=0.0)
        tumor_span_x = tumor_max_x - tumor_min_x
        tumor_span_y = tumor_max_y - tumor_min_y
        tumor_aspect_ratio = tumor_span_x / tumor_span_y if tumor_span_y > 0 else 0.0
        tumor_mean_radius = self._mean(tumor_radial_distances)
        tumor_radius_std = (
            math.sqrt(self._mean([(radius - tumor_mean_radius) ** 2 for radius in tumor_radial_distances]))
            if tumor_radial_distances else 0.0
        )
        all_vessel_positions = [self._cell_com_position(c) for c in all_vascular_cells]
        neo_positions        = [self._cell_com_position(c) for c in neovascular_cells]
    
        if tumor_positions and all_vessel_positions:
            dists_to_any_vessel = [
                self._min_distance_to_points(p, all_vessel_positions)
                for p in tumor_positions
            ]
            mean_dist_to_vessel = self._mean(dists_to_any_vessel)
            min_dist_to_vessel  = min(dists_to_any_vessel)

            prox_thresh = float(getattr(CONFIG, "vessel_proximity_distance", 15.0))
            near_vessel_fraction = sum(
                1 for d in dists_to_any_vessel if d <= prox_thresh
            ) / len(dists_to_any_vessel)
        else:
            mean_dist_to_vessel = min_dist_to_vessel = near_vessel_fraction = 0.0

        if tumor_positions and neo_positions:
            dists_to_neo = [
                self._min_distance_to_points(p, neo_positions)
                for p in tumor_positions
            ]
            mean_dist_to_sprout = self._mean(dists_to_neo)
            min_dist_to_sprout  = min(dists_to_neo)
        else:
            mean_dist_to_sprout = min_dist_to_sprout = 0.0
        proximity_threshold = float(CONFIG.vessel_proximity_distance)
        hypoxic_core_offset = self._euclidean_distance(tumor_centroid, hypoxic_core_centroid)
        hypoxic_mean_radius = self._mean(hypoxic_radial_distances)
        necrotic_mean_radius = self._mean(necrotic_radial_distances)
        normal_radial_distances = self._radial_distances(normal_cells, tumor_centroid)
        normal_mean_radius = self._mean(normal_radial_distances)
        inner_fraction = float(getattr(CONFIG, "hypoxic_core_inner_fraction", 0.45))
        inner_radius   = inner_fraction * tumor_effective_radius
        #inner_core_radius = float(CONFIG.hypoxic_core_inner_fraction) * tumor_effective_radius
        hypoxic_inner_fraction = (
            sum(1 for radius in hypoxic_radial_distances if radius <= inner_radius) / len(hypoxic_radial_distances)
            if hypoxic_radial_distances else 0.0
        )
        necrotic_inner_fraction = (
            sum(1 for radius in necrotic_radial_distances if radius <= inner_radius) / len(necrotic_radial_distances)
            if necrotic_radial_distances else 0.0
        )

        metrics = {
            "mcs": int(mcs),
            "phase": self._timescale_phase_label(mcs),
            "is_growth_window": int(self._is_growth_window(mcs)),
            "tumor_cells": len(tumor_like_cells),
            "normal_cells": len(normal_cells),
            "hypoxic_cells": len(hypoxic_cells),
            "necrotic_cells": len(necrotic_cells),
            "tumor_like_cells": len(tumor_like_cells),
            "hypoxic_fraction": len(hypoxic_cells) / len(tumor_like_cells) if tumor_like_cells else 0.0,
            "necrotic_fraction": len(necrotic_cells) / len(tumor_like_cells) if tumor_like_cells else 0.0,
            "avg_tumor_volume": self._mean(tumor_volumes),
            "avg_normal_volume": self._mean(normal_volumes),
            "avg_hypoxic_volume": self._mean(hypoxic_volumes),
            "avg_necrotic_volume": self._mean(necrotic_volumes),
            "total_tumor_volume": sum(tumor_volumes),
            "avg_tumor_target_volume": self._mean(tumor_target_volumes),
            "tumor_effective_radius": tumor_effective_radius,
            "tumor_mean_radius": tumor_mean_radius,
            "tumor_radius_std": tumor_radius_std,
            "tumor_radial_cv": tumor_radius_std / tumor_mean_radius if tumor_mean_radius > 0 else 0.0,
            "tumor_span_x": tumor_span_x,
            "tumor_span_y": tumor_span_y,
            "tumor_aspect_ratio": tumor_aspect_ratio,
            "mean_tumor_to_vessel_distance": mean_dist_to_vessel,
            "min_tumor_to_vessel_distance": min_dist_to_vessel ,
            "tumor_near_vessel_fraction": near_vessel_fraction,
            "hypoxic_core_offset": hypoxic_core_offset,
            "hypoxic_mean_radius": hypoxic_mean_radius,
            "necrotic_mean_radius": necrotic_mean_radius,
            "hypoxic_inner_fraction": hypoxic_inner_fraction,
            "necrotic_inner_fraction": necrotic_inner_fraction,
            "normal_fraction":    len(normal_cells) / len(tumor_like_cells) if tumor_like_cells else 0.0,
            "vessel_coverage_gap": near_vessel_fraction - (len(normal_cells) / len(tumor_like_cells) if tumor_like_cells else 0.0),
            "normal_mean_radius": normal_mean_radius,
            "mean_dist_to_sprout": mean_dist_to_sprout,
            "min_dist_to_sprout":  min_dist_to_sprout,
        }

        if CONFIG.monitor_include_growth_rates:
            metrics.update({
                "avg_tumor_volume_growth_rate": self._mean_delta_rate(
                    current_tumor_volumes,
                    self._previous_tumor_volumes,
                    delta_mcs,
                ),
                "avg_tumor_target_growth_rate": self._mean_delta_rate(
                    current_tumor_targets,
                    self._previous_tumor_targets,
                    delta_mcs,
                ),
                "total_tumor_volume_growth_rate": self._total_delta_rate(
                    current_tumor_volumes,
                    self._previous_tumor_volumes,
                    delta_mcs,
                ),
            })

        if CONFIG.monitor_include_vascular_metrics:
            metrics.update({
                "endothelial_cells": len(endothelial_cells),
                "neovascular_cells": len(neovascular_cells),
                "vascular_like_cells": len(vascular_like_cells),
                "active_neovascular_cells": len(active_neovascular_cells),
                "vascular_cells": len(vascular_cells),
                "inactive_neovascular_cells": len(inactive_neovascular_cells),
                "avg_active_neovascular_volume": self._mean(active_neovascular_volumes),
                "avg_inactive_neovascular_volume": self._mean(inactive_neovascular_volumes),
                "avg_neovascular_volume": self._mean(neovascular_volumes),
                "avg_endothelial_volume": self._mean(endothelial_volumes),
                "avg_parent_vessel_volume": self._mean(parent_vessel_volumes),
            })
            if CONFIG.monitor_include_growth_rates:
                metrics["avg_endothelial_volume_growth_rate"] = self._mean_delta_rate(
                    current_endothelial_volumes,
                    self._previous_endothelial_volumes,
                    delta_mcs,
                )
                metrics["avg_parent_vessel_volume_growth_rate"] = self._mean_delta_rate(
                    current_vascular_volumes,
                    self._previous_vascular_volumes,
                    delta_mcs,
                )

        if CONFIG.monitor_include_field_means:
            oxygen_field = self.field.Oxygen
            vegf1_field = self.field.VEGF1
            vegf2_field = self.field.VEGF2
            metrics.update({
                "mean_tumor_oxygen": self._mean([self._field_at_com(oxygen_field, cell) for cell in tumor_like_cells]),
                "mean_tumor_vegf": self._mean([self._field_at_com(vegf2_field, cell) for cell in tumor_like_cells]),
                "mean_tumor_vegf1": self._mean([self._field_at_com(vegf1_field, cell) for cell in tumor_like_cells]),
                "mean_tumor_vegf2": self._mean([self._field_at_com(vegf2_field, cell) for cell in tumor_like_cells]),
            })
            if CONFIG.monitor_include_vascular_metrics:
                metrics.update({
                    "mean_endothelial_oxygen": self._mean([
                        self._field_at_com(oxygen_field, cell) for cell in endothelial_cells
                    ]),
                    "mean_endothelial_vegf": self._mean([
                        self._field_at_com(vegf2_field, cell) for cell in endothelial_cells
                    ]),
                    "mean_vascular_oxygen": self._mean([
                        self._field_at_com(oxygen_field, cell) for cell in vascular_cells
                    ]),
                    "mean_vascular_vegf1": self._mean([
                        self._field_at_com(vegf1_field, cell) for cell in vascular_cells
                    ]),
                    "mean_vascular_vegf2": self._mean([
                        self._field_at_com(vegf2_field, cell) for cell in vascular_cells
                    ]),
                })

        if CONFIG.monitor_include_vascular_metrics:
            metrics.update({
                "tumor_oxygen_per_neovascular_cell": (
                    metrics.get("mean_tumor_oxygen", 0.0) / len(neovascular_cells)
                    if neovascular_cells and CONFIG.monitor_include_field_means else 0.0
                ),
                "tumor_oxygen_per_neovascular_volume": (
                    metrics.get("mean_tumor_oxygen", 0.0) / sum(neovascular_volumes)
                    if neovascular_volumes and CONFIG.monitor_include_field_means else 0.0
                ),
                "hypoxia_per_neovascular_cell": (
                    metrics.get("hypoxic_fraction", 0.0) / len(neovascular_cells)
                    if neovascular_cells else 0.0
                ),
            })

        if CONFIG.enable_hif1a_network:
            mean_tumor_hif1a = self._mean([
                float(cell.dict.get("hif1a", CONFIG.hif1a_initial_value))
                for cell in tumor_like_cells
            ])
            mean_tumor_vegf_drive = self._mean([
                float(cell.dict.get("vegf_drive", CONFIG.vegf_drive_basal))
                for cell in tumor_like_cells
            ])
            mean_hif_vegf2_boost = float(CONFIG.hif1a_to_vegf2_weight) * mean_tumor_vegf_drive

            metrics.update({
                "mean_tumor_hif1a": mean_tumor_hif1a,
                "mean_tumor_vegf_drive": mean_tumor_vegf_drive,
                "mean_hif_vegf2_boost": mean_hif_vegf2_boost,
            })

            if CONFIG.monitor_include_field_means and CONFIG.monitor_include_vascular_metrics:
                metrics["mean_endothelial_effective_vegf2"] = (
                    metrics.get("mean_endothelial_vegf", 0.0) + mean_hif_vegf2_boost
                )

        self._previous_sample_mcs = mcs
        self._previous_tumor_volumes = current_tumor_volumes
        self._previous_tumor_targets = current_tumor_targets
        self._previous_endothelial_volumes = current_endothelial_volumes
        self._previous_vascular_volumes = current_vascular_volumes
        return metrics

    def _print_metrics(self, metrics):
        line = (
            f"[Angiogenesis][Monitor][MCS={metrics['mcs']}] "
            f"phase={metrics['phase']} "
            f"tumor_like={metrics['tumor_like_cells']} hypoxic_fraction={metrics['hypoxic_fraction']:.3f} "
            f"necrotic_fraction={metrics['necrotic_fraction']:.3f} avg_tumor_volume={metrics['avg_tumor_volume']:.3f}"
        )

        if CONFIG.monitor_include_growth_rates:
            line += (
                f" avg_dV={metrics['avg_tumor_volume_growth_rate']:.3f}/MCS"
                f" avg_dTargetV={metrics['avg_tumor_target_growth_rate']:.3f}/MCS"
            )

        if CONFIG.monitor_include_field_means:
            line += (
                f" mean_O2={metrics['mean_tumor_oxygen']:.3f}"
                f" mean_VEGF2={metrics['mean_tumor_vegf2']:.3f}"
            )

        if CONFIG.enable_hif1a_network:
            line += (
                f" mean_HIF1a={metrics['mean_tumor_hif1a']:.3f}"
                f" mean_VEGF_drive={metrics['mean_tumor_vegf_drive']:.3f}"
                f" HIF_boost={metrics['mean_hif_vegf2_boost']:.3f}"
            )

        if CONFIG.monitor_include_vascular_metrics:
            line += f" vascular_like={metrics['vascular_like_cells']}"

        print(line)

    def step(self, mcs):
        if not CONFIG.enable_python_monitoring:
            return
        if mcs % CONFIG.monitor_frequency != 0:
            return

        metrics = self._collect_metrics(mcs)
        self.metric_history.append(metrics)

        if CONFIG.monitor_to_csv:
            self._write_metrics_to_csv(metrics)
        if CONFIG.monitor_to_console:
            self._print_metrics(metrics)


class ReporterSteppable(BaseModelSteppable):
    NORMAL: int
    HYPOXIC: int
    NECROTIC: int
    ACTIVENEOVASCULAR: int
    VASCULAR: int
    INACTIVENEOVASCULAR: int

    def step(self, mcs):
        if mcs % CONFIG.report_frequency != 0:
            return

        normal_count = len(list(self.cell_list_by_type(self.NORMAL)))
        hypoxic_count = len(list(self.cell_list_by_type(self.HYPOXIC)))
        necrotic_count = len(list(self.cell_list_by_type(self.NECROTIC)))
        active_count = len(list(self.cell_list_by_type(self.ACTIVENEOVASCULAR)))
        vascular_count = len(list(self.cell_list_by_type(self.VASCULAR)))
        inactive_count = len(list(self.cell_list_by_type(self.INACTIVENEOVASCULAR)))

        print(
            f"[Angiogenesis][MCS={mcs}] phase={self._timescale_phase_label(mcs)} "
            f"normal={normal_count} hypoxic={hypoxic_count} necrotic={necrotic_count} "
            f"active_neo={active_count} vascular={vascular_count} inactive_neo={inactive_count}"
        )
