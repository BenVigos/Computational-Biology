from cc3d.core.PySteppables import *

# Using the standard base but ensuring we use the direct division call
class SimulationSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # 1. ECM
        for x in range(0, 200, 5):
            for y in range(0, 200, 5):
                cell = self.new_cell(self.ECM)
                self.cell_field[x:x+5, y:y+5, 0] = cell
                cell.targetVolume = 25
                cell.lambdaVolume = 4.0

        # 2. Tumor
        cell = self.new_cell(self.TUMOR)
        self.cell_field[95:105, 95:105, 0] = cell
        cell.targetVolume = 100 
        cell.lambdaVolume = 2.0

        # 3. Vessel
        cell = self.new_cell(self.ENDOTHELIAL)
        self.cell_field[5:15, 20:180, 0] = cell
        cell.targetVolume = 1600
        cell.lambdaVolume = 2.0

    def step(self, mcs):
        oxy_field = self.field.Oxygen
        vegf_field = self.field.VEGF
        mmp_field = self.field.MMP

        cells_to_divide = []

        for cell in self.cell_list:
            if not cell: continue
            x, y = int(cell.xCOM), int(cell.yCOM)
            
            if cell.type in [self.TUMOR, self.HYPOXIC]:
                oxy = oxy_field[x, y, 0]
                if oxy > 0.15:
                    cell.targetVolume += 0.5
                    cell.type = self.TUMOR
                else:
                    cell.type = self.HYPOXIC
                    vegf_field[x, y, 0] += 0.2
                if oxy < 0.02:
                    cell.targetVolume -= 1.0

            if cell.type == self.ENDOTHELIAL:
                mmp_field[x, y, 0] += 0.8
                if mcs % 20 == 0:
                    cell.targetVolume += 10

            if cell.type != self.ECM and cell.targetVolume > 60:
                cells_to_divide.append(cell)

        # DIRECT DIVISION CALL
        for cell in cells_to_divide:
            # This is the most primitive division call, bypassing helper attributes
            self.divide_cell_orientation_vector_based(cell, 1, 0, 0)

        # ECM Melting
        for cell in self.cell_list_by_type(self.ECM):
            mmp_val = mmp_field[int(cell.xCOM), int(cell.yCOM), 0]
            v_melt = (2.0 * mmp_val) / (0.5 + mmp_val)
            if v_melt > 0.1:
                cell.targetVolume -= v_melt
                if cell.targetVolume < 1:
                    cell.targetVolume = 0

    def update_attributes(self):
        self.parent_cell.targetVolume /= 2.0
        # This handles the cloning of properties to the new cell
        self.child_cell.targetVolume = self.parent_cell.targetVolume
        self.child_cell.lambdaVolume = self.parent_cell.lambdaVolume
        self.child_cell.type = self.parent_cell.type