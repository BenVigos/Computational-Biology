# High-level steppables for the Angiogenesis simulation.
# This file contains the runtime steppable classes that control initialization
# and dynamics of cells and fields (oxygen, VEGF). Do not change how CC3D
# injects cell type attributes — they come from Angiogenesis.xml (e.g. Medium,
# Tumor, Hypoxic, Endothelial, BloodVessel) and become available as
# self.MEDIUM, self.TUMOR, etc. at runtime.

import csv
import pathlib

from cc3d.core.PySteppables import *

from AngiogenesisConfig import CONFIG


# --- Type hints for editor / linters -------------------------------------------------
# These are type-only annotations so editors know the expected steppable
# attributes. They do NOT create runtime attributes that would conflict with
# CC3D's injection (annotations do not produce class attributes in Python).
class _CellTypeAttrs:
    MEDIUM: int
    TUMOR: int
    HYPOXIC: int
    ENDOTHELIAL: int
    BLOODVESSEL: int


# --- Base steppable -----------------------------------------------------------------
# BaseModelSteppable contains helper utilities used by the concrete steppables
# below (e.g., field accessors, pixel iteration utilities and a helper that
# enables chemotaxis for endothelial cells).
# Important: cell type IDs (self.TUMOR etc.) are provided by CC3D when the
# steppables are registered. To move the blood vessel, edit the CONFIG values
# in AngiogenesisConfig.py (see details below).
class BaseModelSteppable(_CellTypeAttrs, SteppableBasePy):
    # Cell type attributes (e.g. self.MEDIUM, self.TUMOR) are injected by CC3D
    # from Angiogenesis.xml during steppable initialization.

    def __init__(self, frequency=1):
        # frequency controls how often `step` is called (every `frequency` MCS)
        super().__init__(frequency)
        # internal flag to avoid spamming chemotaxis API warnings
        self._chemotaxis_api_warning_shown = False

    def _clamp(self, value, minimum=0.0, maximum=None):
        # Small utility to clamp floating values used for fields and targets.
        if maximum is not None and value > maximum:
            return maximum
        if value < minimum:
            return minimum
        return value

    def _field_at_com(self, field, cell):
        # Return the field value at a cell's center of mass (COM).
        # We round the COM and clamp to lattice bounds before indexing.
        x = min(max(int(round(cell.xCOM)), 0), self.dim.x - 1)
        y = min(max(int(round(cell.yCOM)), 0), self.dim.y - 1)
        z = min(max(int(round(cell.zCOM)), 0), self.dim.z - 1)
        return float(field[x, y, z])

    def _iter_cell_pixels(self, cell):
        # Yield Pixel objects for every pixel belonging to `cell`.
        for pixel_data in self.get_cell_pixel_list(cell):
            yield pixel_data.pixel

    def _add_to_cell_field(self, field, cell, delta):
        # Add `delta` to a field at all pixels occupied by `cell`.
        # The field values are clamped using CONFIG.field_max_value.
        for pixel in self._iter_cell_pixels(cell):
            current_value = float(field[pixel.x, pixel.y, pixel.z])
            field[pixel.x, pixel.y, pixel.z] = self._clamp(
                current_value + delta,
                minimum=0.0,
                maximum=CONFIG.field_max_value,
            )

    def _set_endothelial_chemotaxis(self, cell, lambda_value):
        # Configure chemotaxis for an endothelial cell towards the VEGF field.
        # This helper uses the CompuCell3D chemotaxis API when available.
        try:
            chemotaxis_data = None
            try:
                chemotaxis_data = self.chemotaxisPlugin.getChemotaxisData(cell, "VEGF")
            except Exception:
                chemotaxis_data = None

            if chemotaxis_data is None:
                chemotaxis_data = self.chemotaxisPlugin.addChemotaxisData(cell, "VEGF")

            chemotaxis_data.setLambda(lambda_value)

            # `assignChemotactTowardsVectorTypes` accepts a list of cell type
            # IDs that define where the vector direction points. Medium is used
            # here as the reference type. This call is best-effort and wrapped
            # in try/except because older/newer CC3D APIs differ.
            try:
                chemotaxis_data.assignChemotactTowardsVectorTypes([self.MEDIUM])
            except Exception:
                pass
        except Exception:
            if not self._chemotaxis_api_warning_shown:
                print("[Angiogenesis] Chemotaxis API unavailable; keeping chemotaxis disabled.")
                self._chemotaxis_api_warning_shown = True

    def _mean(self, values):
        # Safe mean helper for metric aggregation.
        return sum(values) / len(values) if values else 0.0


# --- ConstraintInitializerSteppable -------------------------------------------------
# Responsible for creating the initial geometry of the simulation:
# - places a blood vessel region (a rectangular block) defined by
#   CONFIG.vessel_x_min..vessel_x_max and CONFIG.vessel_y_min..vessel_y_max
# - seeds endothelial tip cells at positions defined by CONFIG.tip_cell_y_values
# - seeds tumor cells inside CONFIG tumor rectangle
# To change the blood vessel location, edit the CONFIG values in
# `AngiogenesisConfig.py` (see further below for exact names and an example).
class ConstraintInitializerSteppable(BaseModelSteppable):
    def start(self):
        # Create a blood vessel cell and fill a rectangular region with it.
        vessel = self.new_cell(self.BLOODVESSEL)
        self.cell_field[
            CONFIG.vessel_x_min:CONFIG.vessel_x_max,
            CONFIG.vessel_y_min:CONFIG.vessel_y_max,
            0,
        ] = vessel
        # Give the vessel a target volume consistent with its area.
        vessel.targetVolume = (CONFIG.vessel_x_max - CONFIG.vessel_x_min) * (
            CONFIG.vessel_y_max - CONFIG.vessel_y_min
        )
        vessel.lambdaVolume = CONFIG.blood_vessel_lambda_volume
        vessel.dict["HIF"] = 0.0

        # Place endothelial tip cells at configured Y positions. Each tip cell is
        # a small vertical strip at X between tip_cell_x_min..tip_cell_x_max
        for y_min in CONFIG.tip_cell_y_values:
            endothelial = self.new_cell(self.ENDOTHELIAL)
            self.cell_field[
                CONFIG.tip_cell_x_min:CONFIG.tip_cell_x_max,
                y_min:y_min + CONFIG.tip_cell_height,
                0,
            ] = endothelial
            endothelial.targetVolume = CONFIG.default_target_volume
            endothelial.lambdaVolume = CONFIG.endothelial_lambda_volume
            endothelial.dict["HIF"] = 0.0
            if CONFIG.enable_chemotaxis:
                # Enable chemotaxis towards VEGF (if the model preset enables it)
                self._set_endothelial_chemotaxis(endothelial, CONFIG.chemotaxis_lambda)

        # Seed tumor cells inside the configured tumor rectangle.
        for x_min in range(CONFIG.tumor_x_min, CONFIG.tumor_x_max, CONFIG.tumor_seed_size):
            for y_min in range(CONFIG.tumor_y_min, CONFIG.tumor_y_max, CONFIG.tumor_seed_size):
                tumor = self.new_cell(self.TUMOR)
                self.cell_field[
                    x_min:min(x_min + CONFIG.tumor_seed_size, CONFIG.tumor_x_max),
                    y_min:min(y_min + CONFIG.tumor_seed_size, CONFIG.tumor_y_max),
                    0,
                ] = tumor
                tumor.targetVolume = CONFIG.default_target_volume
                tumor.lambdaVolume = CONFIG.default_lambda_volume
                tumor.dict["HIF"] = 0.0

        print(f"[Angiogenesis] Initialised preset '{CONFIG.preset_name}'.")


# --- FieldDynamicsSteppable ---------------------------------------------------------
# Handles field updates and cell-level state changes each MCS:
# - oxygen supply from vessels and oxygen uptake by tumor cells
# - HIF update and hypoxia switching for tumor cells
# - VEGF secretion by hypoxic tumor cells
class FieldDynamicsSteppable(BaseModelSteppable):
    def step(self, mcs):
        oxygen_field = self.field.Oxygen
        vegf_field = self.field.VEGF

        # Add oxygen to the tissue from blood vessels if enabled
        if CONFIG.enable_oxygen:
            for vessel in self.cell_list_by_type(self.BLOODVESSEL):
                self._add_to_cell_field(oxygen_field, vessel, CONFIG.oxygen_supply_rate)

        # Tumor cells (normal or hypoxic) consume oxygen
        if CONFIG.enable_oxygen_uptake:
            for tumor_cell in self.cell_list_by_type(self.TUMOR, self.HYPOXIC):
                for pixel in self._iter_cell_pixels(tumor_cell):
                    oxygen_value = float(oxygen_field[pixel.x, pixel.y, pixel.z])
                    oxygen_field[pixel.x, pixel.y, pixel.z] = self._clamp(
                        oxygen_value - CONFIG.oxygen_uptake_rate * oxygen_value,
                        minimum=0.0,
                        maximum=CONFIG.field_max_value,
                    )

        # Update HIF and hypoxia state, and possibly secrete VEGF
        for tumor_cell in self.cell_list_by_type(self.TUMOR, self.HYPOXIC):
            local_oxygen = self._field_at_com(oxygen_field, tumor_cell)

            if CONFIG.enable_hif:
                # HIF dynamics (if enabled) follow a simple ODE discretisation
                hif_value = float(tumor_cell.dict.get("HIF", 0.0))
                tumor_cell.dict["HIF"] = self._clamp(
                    hif_value + CONFIG.hif_alpha - CONFIG.hif_beta * local_oxygen * hif_value,
                    minimum=0.0,
                )
            else:
                tumor_cell.dict["HIF"] = 0.0

            # Switch tumor cell phenotype to hypoxic based on oxygen threshold
            if CONFIG.enable_hypoxia_switch:
                tumor_cell.type = self.HYPOXIC if local_oxygen < CONFIG.oxygen_hypoxia_threshold else self.TUMOR

            # Hypoxic tumor cells can secrete VEGF
            if CONFIG.enable_vegf_secretion and tumor_cell.type == self.HYPOXIC:
                secretion_rate = CONFIG.vegf_secretion_rate
                if CONFIG.enable_hif:
                    secretion_rate *= 1.0 + CONFIG.hif_to_vegf_gain * float(tumor_cell.dict.get("HIF", 0.0))
                self._add_to_cell_field(vegf_field, tumor_cell, secretion_rate)


# --- GrowthSteppable ----------------------------------------------------------------
# Updates cell target volumes (growth) for tumor and endothelial cells and
# ensures endothelial chemotaxis is set correctly each MCS.
class GrowthSteppable(BaseModelSteppable):
    def step(self, mcs):
        oxygen_field = self.field.Oxygen
        vegf_field = self.field.VEGF

        # Tumor growth depends on local oxygen availability
        if CONFIG.enable_tumor_growth:
            for tumor_cell in self.cell_list_by_type(self.TUMOR, self.HYPOXIC):
                local_oxygen = self._field_at_com(oxygen_field, tumor_cell)
                growth_rate = CONFIG.tumor_growth_gm * local_oxygen / (
                    CONFIG.oxygen_mm_constant + local_oxygen + 1e-12
                )


                if CONFIG.enable_apoptosis and local_oxygen < CONFIG.oxygen_death_threshold:
                    growth_rate -= CONFIG.apoptosis_rate

                tumor_cell.targetVolume = self._clamp(
                    tumor_cell.targetVolume + growth_rate,
                    minimum=0.0,
                )

        # Endothelial growth is triggered by VEGF
        if CONFIG.enable_endothelial_growth:
            for endothelial_cell in self.cell_list_by_type(self.ENDOTHELIAL):
                local_vegf = self._field_at_com(vegf_field, endothelial_cell)
                if local_vegf < CONFIG.vegf_activation_threshold:
                    continue

                growth_rate = CONFIG.endothelial_growth_gv * local_vegf / (
                    CONFIG.vegf_half_max_scale * CONFIG.vegf_activation_threshold + local_vegf + 1e-12
                )
                endothelial_cell.targetVolume += growth_rate

        # Ensure chemotaxis is configured each MCS if enabled
        chemotaxis_lambda = CONFIG.chemotaxis_lambda if CONFIG.enable_chemotaxis else 0.0
        for endothelial_cell in self.cell_list_by_type(self.ENDOTHELIAL):
            self._set_endothelial_chemotaxis(endothelial_cell, chemotaxis_lambda)


# --- MitosisSteppable ----------------------------------------------------------------
# Handles cell division (mitosis) when enabled in CONFIG. Uses dynamic
# getattr lookups for cell type attributes to avoid static-analysis warnings.
class MitosisSteppable(_CellTypeAttrs, MitosisSteppableBase):
    def __init__(self, frequency=1):
        super().__init__(frequency)

    def step(self, mcs):
        if not CONFIG.enable_mitosis:
            return

        # Use getattr so static linters don't complain; at runtime CC3D will
        # have injected concrete integer IDs for these attributes.
        tumor_type = getattr(self, "TUMOR")
        hypoxic_type = getattr(self, "HYPOXIC")
        endothelial_type = getattr(self, "ENDOTHELIAL")

        cells_to_divide = []
        for cell in self.cell_list_by_type(tumor_type, hypoxic_type, endothelial_type):
            division_threshold = (
                CONFIG.endothelial_division_volume
                if cell.type == endothelial_type
                else CONFIG.tumor_division_volume
            )
            if cell.volume >= division_threshold:
                cells_to_divide.append(cell)

        for cell in cells_to_divide:
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        # Called during mitosis to set attributes for parent & child cells.
        self.parent_cell.targetVolume /= 2.0
        self.clone_parent_2_child()
        self.child_cell.type = self.parent_cell.type


# --- MonitoringSteppable -------------------------------------------------------------
# Samples model-level metrics in Python so you can track how the simulation
# changes over time without hard-coding analysis outside the CC3D run.
# Metrics are stored in `self.metric_history`, and can also be printed and/or
# written to CSV depending on the CONFIG monitoring toggles.
class MonitoringSteppable(BaseModelSteppable):
    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.metric_history = []
        self._previous_sample_mcs = None
        self._previous_tumor_volumes = {}
        self._previous_tumor_targets = {}
        self._previous_endothelial_volumes = {}
        self._csv_path = None

    def start(self):
        if not CONFIG.enable_python_monitoring or not CONFIG.monitor_to_csv:
            return

        output_dir = pathlib.Path(__file__).resolve().parents[1] / CONFIG.monitor_output_dir
        output_dir.mkdir(parents=True, exist_ok=True)
        self._csv_path = output_dir / CONFIG.monitor_output_filename

        # Truncate the file at the beginning of each run so the CSV matches the
        # current simulation configuration.
        self._csv_path.write_text("", encoding="utf-8")
        print(f"[Angiogenesis] Monitoring CSV -> {self._csv_path}")

    def _cell_value_map(self, cells, *, attribute_name):
        value_map = {}
        for cell in cells:
            cell_id = getattr(cell, "id", None)
            if cell_id is None:
                continue
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
        tumor_cells = list(self.cell_list_by_type(self.TUMOR, self.HYPOXIC))
        normoxic_tumor_cells = list(self.cell_list_by_type(self.TUMOR))
        hypoxic_cells = list(self.cell_list_by_type(self.HYPOXIC))
        endothelial_cells = list(self.cell_list_by_type(self.ENDOTHELIAL))

        delta_mcs = 0 if self._previous_sample_mcs is None else mcs - self._previous_sample_mcs
        tumor_volumes = [float(cell.volume) for cell in tumor_cells]
        tumor_target_volumes = [float(cell.targetVolume) for cell in tumor_cells]
        endothelial_volumes = [float(cell.volume) for cell in endothelial_cells]

        current_tumor_volumes = self._cell_value_map(tumor_cells, attribute_name="volume")
        current_tumor_targets = self._cell_value_map(tumor_cells, attribute_name="targetVolume")
        current_endothelial_volumes = self._cell_value_map(endothelial_cells, attribute_name="volume")

        metrics = {
            "mcs": int(mcs),
            "tumor_cells": len(normoxic_tumor_cells),
            "hypoxic_cells": len(hypoxic_cells),
            "tumor_like_cells": len(tumor_cells),
            "hypoxic_fraction": len(hypoxic_cells) / len(tumor_cells) if tumor_cells else 0.0,
            "avg_tumor_volume": self._mean(tumor_volumes),
            "total_tumor_volume": sum(tumor_volumes),
            "avg_tumor_target_volume": self._mean(tumor_target_volumes),
            "mean_hif": self._mean([float(cell.dict.get("HIF", 0.0)) for cell in tumor_cells]),
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

        if CONFIG.monitor_include_endothelial_metrics:
            metrics.update({
                "endothelial_cells": len(endothelial_cells),
                "avg_endothelial_volume": self._mean(endothelial_volumes),
            })
            if CONFIG.monitor_include_growth_rates:
                metrics["avg_endothelial_volume_growth_rate"] = self._mean_delta_rate(
                    current_endothelial_volumes,
                    self._previous_endothelial_volumes,
                    delta_mcs,
                )

        if CONFIG.monitor_include_field_means:
            oxygen_field = self.field.Oxygen
            vegf_field = self.field.VEGF
            metrics.update({
                "mean_tumor_oxygen": self._mean([self._field_at_com(oxygen_field, cell) for cell in tumor_cells]),
                "mean_tumor_vegf": self._mean([self._field_at_com(vegf_field, cell) for cell in tumor_cells]),
            })
            if CONFIG.monitor_include_endothelial_metrics:
                metrics.update({
                    "mean_endothelial_oxygen": self._mean([
                        self._field_at_com(oxygen_field, cell) for cell in endothelial_cells
                    ]),
                    "mean_endothelial_vegf": self._mean([
                        self._field_at_com(vegf_field, cell) for cell in endothelial_cells
                    ]),
                })

        self._previous_sample_mcs = mcs
        self._previous_tumor_volumes = current_tumor_volumes
        self._previous_tumor_targets = current_tumor_targets
        self._previous_endothelial_volumes = current_endothelial_volumes
        return metrics

    def _print_metrics(self, metrics):
        line = (
            f"[Angiogenesis][Monitor][MCS={metrics['mcs']}] "
            f"tumor_like={metrics['tumor_like_cells']} hypoxic_fraction={metrics['hypoxic_fraction']:.3f} "
            f"avg_tumor_volume={metrics['avg_tumor_volume']:.3f}"
        )

        if CONFIG.monitor_include_growth_rates:
            line += (
                f" avg_dV={metrics['avg_tumor_volume_growth_rate']:.3f}/MCS"
                f" avg_dTargetV={metrics['avg_tumor_target_growth_rate']:.3f}/MCS"
            )

        if CONFIG.monitor_include_field_means:
            line += (
                f" mean_O2={metrics['mean_tumor_oxygen']:.3f}"
                f" mean_VEGF={metrics['mean_tumor_vegf']:.3f}"
            )

        if CONFIG.monitor_include_endothelial_metrics:
            line += f" endothelial={metrics['endothelial_cells']}"

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


# --- ReporterSteppable ---------------------------------------------------------------
# Periodically prints simple diagnostics (counts of cell types and mean HIF).
class ReporterSteppable(BaseModelSteppable):
    def step(self, mcs):
        if mcs % CONFIG.report_frequency != 0:
            return

        tumor_count = len(list(self.cell_list_by_type(self.TUMOR)))
        hypoxic_count = len(list(self.cell_list_by_type(self.HYPOXIC)))
        endothelial_count = len(list(self.cell_list_by_type(self.ENDOTHELIAL)))
        mean_hif = 0.0
        tumor_like_cells = list(self.cell_list_by_type(self.TUMOR, self.HYPOXIC))
        if tumor_like_cells:
            mean_hif = sum(float(cell.dict.get("HIF", 0.0)) for cell in tumor_like_cells) / len(tumor_like_cells)

        print(
            f"[Angiogenesis][MCS={mcs}] tumor={tumor_count} hypoxic={hypoxic_count} "
            f"endothelial={endothelial_count} mean_HIF={mean_hif:.3f}"
        )
