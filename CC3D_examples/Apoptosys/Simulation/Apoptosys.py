
from cc3d import CompuCellSetup
        

from ApoptosysSteppables import ApoptosysSteppable

CompuCellSetup.register_steppable(steppable=ApoptosysSteppable(frequency=1))


CompuCellSetup.run()
