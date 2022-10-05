
from cc3d.core.PySteppables import *

# todo: make this robust and non-local
import sys
sys.path.extend(['D:\\modeling\\PhenoCellPy', 'D:/modeling/PhenoCellPy'])

import Phenotypes as pheno


class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        
        side = 6
        
        x = self.dim.x//2 - side//2
        y = self.dim.x//2 - side//2
        
        cell = self.new_cell(self.CELL)
        self.cell_field[x:x + side - 1, y:y + side - 1, 0] = cell

        dt = 1  # 1 min/mcs

        ki67_basic = pheno.cycles.Ki67Basic(dt=dt, target_volumes=[side*side, side*side],
                                            volumes=[side*side, side*side])

        for cell in self.cell_list:
            cell.targetVolume = side*side
            cell.lambdaVolume = 2.0
            # print(hasattr(cell, "__dict__"))
            pheno.utils.add_cycle_to_object(cell, ki67_basic)
            # print(hasattr(cell, "cycle"))

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):
        
        if not mcs:
            side = 6
            dt = 1  # 1 min/mcs
            ki67_basic = pheno.cycles.Ki67Basic(dt=1, target_volumes=[side*side, side*side],
                                            volumes=[side*side, side*side])

            for cell in self.cell_list:
                # pheno.utils.add_cycle_to_object(cell, ki67_basic)
                # setattr(cell, "cycle", ki67_basic)
                cell.cycle = ki67_basic
                print(hasattr(cell, "cycle"))

        cells_to_divide = []
        for cell in self.cell_list:
            print(hasattr(cell, "cycle"))
            cell.cycle.current_phase.volume = cell.volume
            changed_phase, died, divides = cell.cycle.time_step_cycle()
            cell.targetVolume = cell.cycle.current_phase.volume
            if divides:
                cells_to_divide.append(cell)

        for cell in cells_to_divide:

            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume = 36  # todo: use parameter

        self.clone_parent_2_child()            

        self.parent_cell.cycle.current_phase.target_volume = self.parent_cell.targetVolume
        self.parent_cell.cycle.current_phase.volume = self.parent_cell.targetVolume

        self.child_cell.current_phase.target_volume = self.parent_cell.targetVolume
        self.child_cell.cycle.current_phase.volume = self.parent_cell.targetVolume
