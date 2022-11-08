
from cc3d import CompuCellSetup
        

from NecrosisSteppables import NecrosisSteppable

CompuCellSetup.register_steppable(steppable=NecrosisSteppable(frequency=1))


CompuCellSetup.run()
