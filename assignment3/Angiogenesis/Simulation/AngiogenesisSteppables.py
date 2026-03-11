# Paper-aligned 2D angiogenesis steppables.
# The mechanics mirror `Paper files/Simulation/TumorVasc3DSteppables.py`
# while keeping the lattice 2D and the initial geometry programmatic.

import csv
import pathlib

from cc3d.core.PySteppables import *

from AngiogenesisConfig import CONFIG


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
            cell.targetSurface = CONFIG.vascular_target_surface
            cell.lambdaSurface = CONFIG.vascular_lambda_surface
        else:
            cell.targetVolume = CONFIG.tumor_target_volume
            cell.lambdaVolume = CONFIG.tumor_lambda_volume
            cell.targetSurface = CONFIG.tumor_target_surface
            cell.lambdaSurface = CONFIG.tumor_lambda_surface


class ConstraintInitializerSteppable(BaseModelSteppable):
    def start(self):
        vessel = self.new_cell(self.VASCULAR)
        self.cell_field[
            CONFIG.vessel_x_min:CONFIG.vessel_x_max,
            CONFIG.vessel_y_min:CONFIG.vessel_y_max,
            0,
        ] = vessel
        self._apply_paper_mechanics(vessel)

        for y_min in CONFIG.tip_cell_y_values:
            tip = self.new_cell(self.INACTIVENEOVASCULAR)
            self.cell_field[
                CONFIG.tip_cell_x_min:CONFIG.tip_cell_x_max,
                y_min:y_min + CONFIG.tip_cell_height,
                0,
            ] = tip
            self._apply_paper_mechanics(tip)

        for x_min in range(CONFIG.tumor_x_min, CONFIG.tumor_x_max, CONFIG.tumor_seed_size):
            for y_min in range(CONFIG.tumor_y_min, CONFIG.tumor_y_max, CONFIG.tumor_seed_size):
                tumor = self.new_cell(self.NORMAL)
                self.cell_field[
                    x_min:min(x_min + CONFIG.tumor_seed_size, CONFIG.tumor_x_max),
                    y_min:min(y_min + CONFIG.tumor_seed_size, CONFIG.tumor_y_max),
                    0,
                ] = tumor
                self._apply_paper_mechanics(tumor)

        print("[Angiogenesis] Initialised paper-aligned 2D geometry.")


class GrowthSteppable(BaseModelSteppable):
    def start(self):
        for cell in self.cellList:
            self._apply_paper_mechanics(cell)

    def _vascular_neighbor_area(self, cell):
        total_area = 0.0
        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
            if neighbor and neighbor.type in self._vascular_type_ids():
                total_area += float(common_surface_area)
        return total_area

    def _vascular_growth_step(self, cell, local_vegf2, neighbor_area_limit):
        if local_vegf2 <= CONFIG.vascular_vegf_activation_threshold:
            return

        if self._vascular_neighbor_area(cell) >= neighbor_area_limit:
            return

        cell.targetVolume += (
            CONFIG.vascular_growth_volume_rate
            * local_vegf2
            / (CONFIG.vascular_growth_denominator + local_vegf2)
        )
        cell.targetSurface += (
            CONFIG.vascular_growth_surface_rate
            * local_vegf2
            / (CONFIG.vascular_growth_denominator + local_vegf2)
        )

    def _tumor_growth_step(self, cell, local_oxygen, mcs):
        if not CONFIG.enable_type_switching:
            return

        if local_oxygen < CONFIG.nutrient_thresh and mcs > CONFIG.tumor_growth_start_mcs:
            cell.type = self.HYPOXIC

        if local_oxygen < CONFIG.necrotic_thresh and mcs > CONFIG.tumor_growth_start_mcs:
            cell.type = self.NECROTIC

        if CONFIG.enable_tumor_growth and mcs > CONFIG.tumor_growth_start_mcs and cell.type in self._growing_tumor_type_ids():
            cell.targetVolume += (
                CONFIG.tumor_growth_volume_rate
                * local_oxygen
                / (CONFIG.tumor_growth_denominator + local_oxygen)
            )
            cell.targetSurface += (
                CONFIG.tumor_growth_surface_rate
                * local_oxygen
                / (CONFIG.tumor_growth_denominator + local_oxygen)
            )

    def _hypoxic_adjustments(self, cell, local_oxygen, mcs):
        if not CONFIG.enable_type_switching:
            return

        if local_oxygen < CONFIG.necrotic_thresh and mcs > CONFIG.tumor_growth_start_mcs:
            cell.type = self.NECROTIC
        elif local_oxygen > CONFIG.nutrient_thresh:
            cell.type = self.NORMAL

    def _necrotic_adjustments(self, cell):
        cell.targetVolume = self._clamp(
            float(cell.targetVolume) - CONFIG.necrotic_volume_loss_rate,
            minimum=0.0,
        )
        cell.lambdaSurface = 0.0

    def step(self, mcs):
        field_neo_vasc = self.field.VEGF2
        field_malig = self.field.Oxygen

        for cell in self.cellList:
            if cell.type == self.INACTIVENEOVASCULAR and CONFIG.enable_vascular_growth:
                self._vascular_growth_step(
                    cell,
                    self._field_at_com(field_neo_vasc, cell),
                    CONFIG.inactive_neighbor_area_limit,
                )

            if cell.type == self.ACTIVENEOVASCULAR and CONFIG.enable_vascular_growth:
                self._vascular_growth_step(
                    cell,
                    self._field_at_com(field_neo_vasc, cell),
                    CONFIG.active_neighbor_area_limit,
                )

            if cell.type in self._growing_tumor_type_ids():
                self._tumor_growth_step(cell, self._field_at_com(field_malig, cell), mcs)

            if cell.type == self.HYPOXIC:
                self._hypoxic_adjustments(cell, self._field_at_com(field_malig, cell), mcs)

            if cell.type == self.NECROTIC:
                self._necrotic_adjustments(cell)


class MitosisSteppable(MitosisSteppableBase):
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
            child_cell.type = self.NORMAL
            child_cell.targetVolume = CONFIG.tumor_target_volume
            child_cell.lambdaVolume = CONFIG.tumor_lambda_volume
            child_cell.targetSurface = CONFIG.tumor_target_surface
            child_cell.lambdaSurface = CONFIG.tumor_lambda_surface
            parent_cell.targetVolume = CONFIG.tumor_target_volume
            parent_cell.lambdaVolume = CONFIG.tumor_lambda_volume
            parent_cell.targetSurface = CONFIG.tumor_target_surface
            parent_cell.lambdaSurface = CONFIG.tumor_lambda_surface

        if parent_cell.type in (self.INACTIVENEOVASCULAR, self.ACTIVENEOVASCULAR):
            child_cell.type = self.ACTIVENEOVASCULAR
            child_cell.targetVolume = CONFIG.vascular_target_volume
            child_cell.lambdaVolume = CONFIG.vascular_lambda_volume
            child_cell.targetSurface = CONFIG.vascular_target_surface
            child_cell.lambdaSurface = CONFIG.vascular_lambda_surface
            parent_cell.targetVolume = CONFIG.vascular_target_volume
            parent_cell.lambdaVolume = CONFIG.vascular_lambda_volume
            parent_cell.targetSurface = CONFIG.vascular_target_surface
            parent_cell.lambdaSurface = CONFIG.vascular_lambda_surface


class MonitoringSteppable(BaseModelSteppable):
    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.metric_history = []
        self._previous_sample_mcs = None
        self._previous_tumor_volumes = {}
        self._previous_tumor_targets = {}
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

        shared_ids = current_values.keys() & previous_values.keys()
        if not shared_ids:
            return 0.0

        return self._mean([
            (current_values[cell_id] - previous_values[cell_id]) / delta_mcs
            for cell_id in shared_ids
        ])

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

        delta_mcs = 0 if self._previous_sample_mcs is None else mcs - self._previous_sample_mcs
        tumor_volumes = [float(cell.volume) for cell in tumor_like_cells]
        tumor_target_volumes = [float(cell.targetVolume) for cell in tumor_like_cells]
        vascular_volumes = [float(cell.volume) for cell in vascular_like_cells]

        current_tumor_volumes = self._cell_value_map(tumor_like_cells, attribute_name="volume")
        current_tumor_targets = self._cell_value_map(tumor_like_cells, attribute_name="targetVolume")
        current_vascular_volumes = self._cell_value_map(vascular_like_cells, attribute_name="volume")

        metrics = {
            "mcs": int(mcs),
            "tumor_cells": len(normal_cells),
            "normal_cells": len(normal_cells),
            "hypoxic_cells": len(hypoxic_cells),
            "necrotic_cells": len(necrotic_cells),
            "tumor_like_cells": len(tumor_like_cells),
            "hypoxic_fraction": len(hypoxic_cells) / len(tumor_like_cells) if tumor_like_cells else 0.0,
            "necrotic_fraction": len(necrotic_cells) / len(tumor_like_cells) if tumor_like_cells else 0.0,
            "avg_tumor_volume": self._mean(tumor_volumes),
            "total_tumor_volume": sum(tumor_volumes),
            "avg_tumor_target_volume": self._mean(tumor_target_volumes),
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
                "endothelial_cells": len(vascular_like_cells),
                "vascular_like_cells": len(vascular_like_cells),
                "active_neovascular_cells": len(active_neovascular_cells),
                "vascular_cells": len(vascular_cells),
                "inactive_neovascular_cells": len(inactive_neovascular_cells),
                "avg_endothelial_volume": self._mean(vascular_volumes),
                "avg_vascular_volume": self._mean(vascular_volumes),
            })
            if CONFIG.monitor_include_growth_rates:
                metrics["avg_endothelial_volume_growth_rate"] = self._mean_delta_rate(
                    current_vascular_volumes,
                    self._previous_vascular_volumes,
                    delta_mcs,
                )
                metrics["avg_vascular_volume_growth_rate"] = metrics["avg_endothelial_volume_growth_rate"]

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
                        self._field_at_com(oxygen_field, cell) for cell in vascular_like_cells
                    ]),
                    "mean_endothelial_vegf": self._mean([
                        self._field_at_com(vegf2_field, cell) for cell in vascular_like_cells
                    ]),
                    "mean_vascular_oxygen": self._mean([
                        self._field_at_com(oxygen_field, cell) for cell in vascular_like_cells
                    ]),
                    "mean_vascular_vegf1": self._mean([
                        self._field_at_com(vegf1_field, cell) for cell in vascular_like_cells
                    ]),
                    "mean_vascular_vegf2": self._mean([
                        self._field_at_com(vegf2_field, cell) for cell in vascular_like_cells
                    ]),
                })

        self._previous_sample_mcs = mcs
        self._previous_tumor_volumes = current_tumor_volumes
        self._previous_tumor_targets = current_tumor_targets
        self._previous_vascular_volumes = current_vascular_volumes
        return metrics

    def _print_metrics(self, metrics):
        line = (
            f"[Angiogenesis][Monitor][MCS={metrics['mcs']}] "
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
            f"[Angiogenesis][MCS={mcs}] normal={normal_count} hypoxic={hypoxic_count} necrotic={necrotic_count} "
            f"active_neo={active_count} vascular={vascular_count} inactive_neo={inactive_count}"
        )
