
from cc3d import CompuCellSetup
        


from Ki67CycleSteppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




# from Ki67CycleSteppables import GrowthSteppable

# CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))




from Ki67CycleSteppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))


CompuCellSetup.run()
