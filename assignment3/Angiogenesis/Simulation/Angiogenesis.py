from cc3d import CompuCellSetup

from AngiogenesisSteppables import (
    ConstraintInitializerSteppable,
    GrowthSteppable,
    MitosisSteppable,
    MonitoringSteppable,
    ReporterSteppable,
)

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=MonitoringSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=ReporterSteppable(frequency=1))

CompuCellSetup.run()
