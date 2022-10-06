from cc3d.cpp.PlayerPython import *
from cc3d import CompuCellSetup

from cc3d.core.PySteppables import *

from numpy import mean

# todo: make this robust and non-local
import sys

# sys.path.extend(['D:\\modeling\\PhenoCellPy', 'D:/modeling/PhenoCellPy'])
sys.path.extend(['D:\\github\\PhenoCellPy', 'D:/github/PhenoCellPy'])

import Phenotypes as pheno


class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.track_cell_level_scalar_attribute(field_name='phase_index_plus_1', attribute_name='phase_index_plus_1')

    def start(self):
        side = 10

        x = self.dim.x // 2 - side // 2
        y = self.dim.x // 2 - side // 2

        cell = self.new_cell(self.CELL)
        self.cell_field[x:x + side, y:y + side, 0] = cell

        dt = 5  # 5 min/mcs

        ki67_basic = pheno.cycles.Ki67Basic(dt=dt, target_volumes=[side * side, side * side],
                                            volumes=[side * side, side * side])

        for cell in self.cell_list:
            cell.targetVolume = side * side
            cell.lambdaVolume = 2.0
            # print(hasattr(cell, "__dict__"))
            pheno.utils.add_cycle_to_CC3D_cell(cell, ki67_basic)
            # print(hasattr(cell, "cycle"))
            cell.dict["phase_index_plus_1"] = cell.dict["cycle"].current_phase.index + 1


class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)
        self.previous_number_cells = 0

        self.plot = True

    def start(self):

        if self.plot:
            self.plot_win_vol = self.add_new_plot_window(title='Average Volume',
                                                         x_axis_title='MonteCarlo Step (MCS)',
                                                         y_axis_title='Variables', x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=True)

            self.plot_win_vol.add_plot("MVol", style='Lines', color='red', size=5)

            self.plot_win_number = self.add_new_plot_window(title='Number of cells',
                                                            x_axis_title='MonteCarlo Step (MCS)',
                                                            y_axis_title='Variables', x_scale_type='linear',
                                                            y_scale_type='linear',
                                                            grid=True)

            self.plot_win_number.add_plot("N", style='Lines', color='red', size=5)

            # initialize setting for Histogram
            self.plot_win_phase = self.add_new_plot_window(title='Number of cells in each phase',
                                                           x_axis_title='MonteCarlo Step (MCS)',
                                                           y_axis_title='Variables', x_scale_type='linear',
                                                           y_scale_type='linear',
                                                           grid=True)

            self.plot_win_phase.add_plot("N", style='Points', color='red', size=20)

    def step(self, mcs):

        if not mcs and self.plot:
            self.plot_win_number.add_data_point("N", mcs, len(self.cell_list))
            self.previous_number_cells = len(self.cell_list)

        elif not mcs % 50 and len(self.cell_list) - self.previous_number_cells > 0 and self.plot:
            self.plot_win_number.add_data_point("N", mcs, len(self.cell_list))
            self.previous_number_cells = len(self.cell_list)

        cells_to_divide = []

        n_zero = 0
        n_one = 0

        volumes = []

        for cell in self.cell_list:
            volumes.append(cell.volume)
            if cell.dict["cycle"].current_phase.index == 0:
                n_zero += 1
            elif cell.dict["cycle"].current_phase.index == 1:
                n_one += 1
            changed_phase, died, divides = cell.dict["cycle"].time_step_cycle()

            if cell.targetVolume < cell.dict["cycle"].current_phase.volume:
                cell.targetVolume = cell.dict["cycle"].current_phase.volume

            if changed_phase:
                cell.dict["phase_index_plus_1"] = cell.dict["cycle"].current_phase.index + 1
                if len(self.cell_list) < 10:
                    print("@@@\nPHASE CHANGE\n@@@")

            if divides:
                if len(self.cell_list) < 10:
                    print(f"@@@\nCELL DIVISION\n@@@\ncell volume={cell.volume}")
                cells_to_divide.append(cell)

        if self.plot:
            self.plot_win_phase.erase_all_data()
            # arguments are (name of the data series, x, y)
            self.plot_win_phase.add_data_point("N", 0, n_zero)
            self.plot_win_phase.add_data_point("N", 1, n_one)

            if not mcs % 50:
                self.plot_win_vol.add_data_point("MVol", mcs, mean(volumes))

        for cell in cells_to_divide:
            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume = 100  # todo: use parameter

        self.clone_parent_2_child()

        self.parent_cell.dict["cycle"].current_phase.target_volume = self.parent_cell.targetVolume
        self.parent_cell.dict["cycle"].current_phase.volume = self.parent_cell.targetVolume
        self.parent_cell.dict["phase_index_plus_1"] = self.parent_cell.dict["cycle"].current_phase.index + 1

        self.child_cell.dict["cycle"].current_phase.target_volume = self.parent_cell.targetVolume
        self.child_cell.dict["cycle"].current_phase.volume = self.parent_cell.targetVolume
        self.child_cell.dict["phase_index_plus_1"] = self.child_cell.dict["cycle"].current_phase.index + 1
        if len(self.cell_list) < 10:
            print("@@@\nCHILD ATTRIBS\n@@@\n", self.child_cell.volume, self.child_cell.dict["cycle"].time_in_cycle,
                  self.child_cell.dict["cycle"].current_phase,
                  self.child_cell.dict["cycle"].current_phase.time_in_phase)
        # self.child_cell.dict["cycle"].time_in_cycle = 0
