from cc3d import CompuCellSetup
from all_values import relaxationmcs
from endothelialSteppables import endothelialSteppable,CellDifferentiationSteppable
from SecretionStep import SecretionSteppables



CompuCellSetup.register_steppable(steppable=endothelialSteppable(frequency=int(relaxationmcs)))
CompuCellSetup.register_steppable(steppable=SecretionSteppables(frequency=int(relaxationmcs //2  )))

CompuCellSetup.register_steppable(steppable=CellDifferentiationSteppable(frequency=1, transition_mcs=10))




CompuCellSetup.run()
