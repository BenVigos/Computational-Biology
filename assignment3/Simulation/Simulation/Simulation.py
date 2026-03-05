from cc3d import CompuCellSetup
from SimulationSteppables import SimulationSteppable

CompuCellSetup.register_steppable(steppable=SimulationSteppable(frequency=1))
CompuCellSetup.run()