from cc3d.core.PySteppables import *

class TumorLogicSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        # Access the Oxygen field
        oxy_field = self.field.Oxygen
        
        for cell in self.cell_list:
            # 1. METABOLISM: Get oxygen at cell center
            local_oxy = oxy_field[cell.x_COM, cell.y_COM, 0]
            
            # 2. GENE REGULATION / TYPE SWITCHING
            # Simple "Hypoxic Switch": If oxygen is low, change color/type
            if local_oxy < 0.2 and cell.type == self.TUMOR:
                cell.type = self.HYPOXIC
            elif local_oxy >= 0.2 and cell.type == self.HYPOXIC:
                cell.type = self.TUMOR
            
            # 3. GROWTH RULE
            # If healthy and fed, increase target volume (prepare to divide)
            if cell.type == self.TUMOR:
                cell.targetVolume += 0.1 
            
            # 4. MITOSIS (Division)
            if cell.targetVolume > 50:
                self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        parent_cell = self.mitosisSteppable.parentCell
        child_cell = self.mitosisSteppable.childCell
        
        # Reset volumes after division
        parent_cell.targetVolume = 25
        child_cell.targetVolume = 25