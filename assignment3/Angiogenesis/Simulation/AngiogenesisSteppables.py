from cc3d.core.PySteppables import *

from AngiogenesisConfig import CONFIG


class _CellTypeAttrs:
    MEDIUM: int
    TUMOR: int
    HYPOXIC: int
    ENDOTHELIAL: int
    BLOODVESSEL: int


class BaseModelSteppable(_CellTypeAttrs, SteppableBasePy):
    # Cell type attributes (e.g. self.MEDIUM, self.TUMOR) are injected by CC3D
    # from Angiogenesis.xml during steppable initialization.

    def __init__(self, frequency=1):
        super().__init__(frequency)
        self._chemotaxis_api_warning_shown = False

    def _clamp(self, value, minimum=0.0, maximum=None):
        if maximum is not None and value > maximum:
            return maximum
        if value < minimum:
            return minimum
        return value

    def _field_at_com(self, field, cell):
        x = min(max(int(round(cell.xCOM)), 0), self.dim.x - 1)
        y = min(max(int(round(cell.yCOM)), 0), self.dim.y - 1)
        z = min(max(int(round(cell.zCOM)), 0), self.dim.z - 1)
        return float(field[x, y, z])

    def _iter_cell_pixels(self, cell):
        for pixel_data in self.get_cell_pixel_list(cell):
            yield pixel_data.pixel

    def _add_to_cell_field(self, field, cell, delta):
        for pixel in self._iter_cell_pixels(cell):
            current_value = float(field[pixel.x, pixel.y, pixel.z])
            field[pixel.x, pixel.y, pixel.z] = self._clamp(
                current_value + delta,
                minimum=0.0,
                maximum=CONFIG.field_max_value,
            )

    def _set_endothelial_chemotaxis(self, cell, lambda_value):
        try:
            chemotaxis_data = None
            try:
                chemotaxis_data = self.chemotaxisPlugin.getChemotaxisData(cell, "VEGF")
            except Exception:
                chemotaxis_data = None

            if chemotaxis_data is None:
                chemotaxis_data = self.chemotaxisPlugin.addChemotaxisData(cell, "VEGF")

            chemotaxis_data.setLambda(lambda_value)

            try:
                chemotaxis_data.assignChemotactTowardsVectorTypes([self.MEDIUM])
            except Exception:
                pass
        except Exception:
            if not self._chemotaxis_api_warning_shown:
                print("[Angiogenesis] Chemotaxis API unavailable; keeping chemotaxis disabled.")
                self._chemotaxis_api_warning_shown = True


class ConstraintInitializerSteppable(BaseModelSteppable):
    def start(self):
        vessel = self.new_cell(self.BLOODVESSEL)
        self.cell_field[
            CONFIG.vessel_x_min:CONFIG.vessel_x_max,
            CONFIG.vessel_y_min:CONFIG.vessel_y_max,
            0,
        ] = vessel
        vessel.targetVolume = (CONFIG.vessel_x_max - CONFIG.vessel_x_min) * (
            CONFIG.vessel_y_max - CONFIG.vessel_y_min
        )
        vessel.lambdaVolume = CONFIG.blood_vessel_lambda_volume
        vessel.dict["HIF"] = 0.0

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
                self._set_endothelial_chemotaxis(endothelial, CONFIG.chemotaxis_lambda)

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


class FieldDynamicsSteppable(BaseModelSteppable):
    def step(self, mcs):
        oxygen_field = self.field.Oxygen
        vegf_field = self.field.VEGF

        if CONFIG.enable_oxygen:
            for vessel in self.cell_list_by_type(self.BLOODVESSEL):
                self._add_to_cell_field(oxygen_field, vessel, CONFIG.oxygen_supply_rate)

        if CONFIG.enable_oxygen_uptake:
            for tumor_cell in self.cell_list_by_type(self.TUMOR, self.HYPOXIC):
                for pixel in self._iter_cell_pixels(tumor_cell):
                    oxygen_value = float(oxygen_field[pixel.x, pixel.y, pixel.z])
                    oxygen_field[pixel.x, pixel.y, pixel.z] = self._clamp(
                        oxygen_value - CONFIG.oxygen_uptake_rate * oxygen_value,
                        minimum=0.0,
                        maximum=CONFIG.field_max_value,
                    )

        for tumor_cell in self.cell_list_by_type(self.TUMOR, self.HYPOXIC):
            local_oxygen = self._field_at_com(oxygen_field, tumor_cell)

            if CONFIG.enable_hif:
                hif_value = float(tumor_cell.dict.get("HIF", 0.0))
                tumor_cell.dict["HIF"] = self._clamp(
                    hif_value + CONFIG.hif_alpha - CONFIG.hif_beta * local_oxygen * hif_value,
                    minimum=0.0,
                )
            else:
                tumor_cell.dict["HIF"] = 0.0

            if CONFIG.enable_hypoxia_switch:
                tumor_cell.type = self.HYPOXIC if local_oxygen < CONFIG.oxygen_hypoxia_threshold else self.TUMOR

            if CONFIG.enable_vegf_secretion and tumor_cell.type == self.HYPOXIC:
                secretion_rate = CONFIG.vegf_secretion_rate
                if CONFIG.enable_hif:
                    secretion_rate *= 1.0 + CONFIG.hif_to_vegf_gain * float(tumor_cell.dict.get("HIF", 0.0))
                self._add_to_cell_field(vegf_field, tumor_cell, secretion_rate)


class GrowthSteppable(BaseModelSteppable):
    def step(self, mcs):
        oxygen_field = self.field.Oxygen
        vegf_field = self.field.VEGF

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

        if CONFIG.enable_endothelial_growth:
            for endothelial_cell in self.cell_list_by_type(self.ENDOTHELIAL):
                local_vegf = self._field_at_com(vegf_field, endothelial_cell)
                if local_vegf < CONFIG.vegf_activation_threshold:
                    continue

                growth_rate = CONFIG.endothelial_growth_gv * local_vegf / (
                    CONFIG.vegf_half_max_scale * CONFIG.vegf_activation_threshold + local_vegf + 1e-12
                )
                endothelial_cell.targetVolume += growth_rate

        chemotaxis_lambda = CONFIG.chemotaxis_lambda if CONFIG.enable_chemotaxis else 0.0
        for endothelial_cell in self.cell_list_by_type(self.ENDOTHELIAL):
            self._set_endothelial_chemotaxis(endothelial_cell, chemotaxis_lambda)


class MitosisSteppable(_CellTypeAttrs, MitosisSteppableBase):
    def __init__(self, frequency=1):
        super().__init__(frequency)

    def step(self, mcs):
        if not CONFIG.enable_mitosis:
            return

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
        self.parent_cell.targetVolume /= 2.0
        self.clone_parent_2_child()
        self.child_cell.type = self.parent_cell.type


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
