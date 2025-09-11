from cc3d import CompuCellSetup
from all_values import relaxationmcs
from CellSteppable import CellSteppable,CellDifferentiationSteppable
from SecretionStep import SecretionSteppables



CompuCellSetup.register_steppable(steppable=CellSteppable(frequency=int(relaxationmcs)))
CompuCellSetup.register_steppable(steppable=SecretionSteppables(frequency=int(relaxationmcs //2  )))

CompuCellSetup.register_steppable(steppable=CellDifferentiationSteppable(frequency=1, transition_mcs=10))




CompuCellSetup.run()
