"""
BSD 3-Clause License

Copyright (c) 2022, Juliano Ferrari Gianlupi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from cc3d.cpp.PlayerPython import *
from cc3d import CompuCellSetup

from cc3d.core.PySteppables import *

from numpy import median, quantile, nan

import sys
from os.path import abspath, dirname, join

# sys.path.extend([abspath("../../..")])  # todo: make this more refined

sys.path.extend(['D:\\modeling\\PhenoCellPy', 'D:/modeling/PhenoCellPy'])
import Phenotypes as pheno


def Ki67pos_transition(*args):
    # print(len(args), print(args))
    # args = [cc3d cell volume, phase's target volume, time in phase, phase duration
    return args[0] >= args[1] and args[2] > args[3]


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

        ki67_basic_modified_transition = pheno.phenotypes.Ki67Basic(dt=dt, target_fluid_fraction=[1, 1],
                                                                    # as the simulated cell "doesn't have" a nucleus
                                                                    # we don't need to give it a volume
                                                                    nuclear_fluid=[0, 0], nuclear_solid=[0, 0],
                                                                    cytoplasm_fluid=[side * side, side * side],
                                                                    cytoplasm_solid=[0, 0],
                                                                    cytoplasm_solid_target=[0, 0],
                                                                    target_cytoplasm_to_nuclear_ratio=[0, 0],
                                                                    transitions_to_next_phase=[None,
                                                                                               Ki67pos_transition],
                                                                    transitions_to_next_phase_args=[None,
                                                                                                    [-9, 1, -9, 1]])

        for cell in self.cell_list:
            cell.targetVolume = side * side
            cell.lambdaVolume = 2.0
            # print(hasattr(cell, "__dict__"))
            pheno.utils.add_phenotype_to_CC3D_cell(cell, ki67_basic_modified_transition)
            # print(hasattr(cell, "phenotype"))
            cell.dict["phase_index_plus_1"] = cell.dict["phenotype"].current_phase.index + 1


class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)
        self.previous_number_cells = 0

        self.plot = True
        self.save = False

        if self.save:
            self.save_loc = dirname(abspath(__file__))

            self.volume_file = open(join(self.save_loc, "volume.dat"), "w+")
            self.volume_file.write("MCS, Median, Min, Max, 10th, 90th, 25th, 75th\n")

            self.time_minus_file = open(join(self.save_loc, "time_in_Ki67-.dat"), "w+")
            self.time_minus_file.write("MCS, Median, Min, Max, 10th, 90th, 25th, 75th\n")

            self.time_plus_file = open(join(self.save_loc, "time_in_Ki67+.dat"), "w+")
            self.time_plus_file.write("MCS, Median, Min, Max, 10th, 90th, 25th, 75th\n")

            self.number_cells_file = open(join(self.save_loc, "number_cells.dat"), "w+")
            self.number_cells_file.write("MCS, N, N+, N-\n")

    def start(self):

        if self.plot:
            self.plot_win_vol = self.add_new_plot_window(title='Volume metrics',
                                                         x_axis_title='MonteCarlo Step (MCS)',
                                                         y_axis_title='Variables', x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=True,
                                                         config_options={"legend": True})

            self.plot_win_vol.add_plot("Median Vol", style='Lines', color='yellow', size=5)
            self.plot_win_vol.add_plot("Minimum Vol", style='Lines', color='blue', size=5)
            self.plot_win_vol.add_plot("Maximum Vol", style='Lines', color='red', size=5)

            self.plot_win_time = self.add_new_plot_window(title='Time spent in each phase',
                                                          x_axis_title='MonteCarlo Step (MCS)',
                                                          y_axis_title='Variables', x_scale_type='linear',
                                                          y_scale_type='linear',
                                                          grid=True,
                                                          config_options={"legend": True})

            self.plot_win_time.add_plot("Median Time in Ki67-", style='Lines', color='yellow', size=5)
            self.plot_win_time.add_plot("Median Time in Ki67+", style='Lines', color='blue', size=5)
            self.plot_win_time.add_plot("Maximum Time in Ki67+", style='Lines', color='red', size=5)

            self.plot_win_number = self.add_new_plot_window(title='Number of cells',
                                                            x_axis_title='MonteCarlo Step (MCS)',
                                                            y_axis_title='Variables', x_scale_type='linear',
                                                            y_scale_type='linear',
                                                            grid=True)
            self.plot_win_number.add_plot("N", style='Lines', color='red', size=5)

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

        time_spent_in_0 = []
        time_spent_in_1 = []

        for cell in self.cell_list:
            # if cell.volume<=90:
            #     print(cell.volume, mcs)
            volumes.append(cell.volume)
            cell.dict["phenotype"].current_phase.simulated_cell_volume = cell.volume
            # print(len(args))
            if cell.dict["phenotype"].current_phase.index == 0:
                n_zero += 1
                time_spent_in_0.append(cell.dict["phenotype"].current_phase.time_in_phase)
            elif cell.dict["phenotype"].current_phase.index == 1:
                n_one += 1
                time_spent_in_1.append(cell.dict["phenotype"].current_phase.time_in_phase)
                # args = [cc3d cell volume, phase's target volume, time in phase, phase duration
                args = [
                    cell.volume,
                    cell.targetVolume,
                    cell.dict["phenotype"].current_phase.time_in_phase + cell.dict["phenotype"].dt,
                    cell.dict["phenotype"].current_phase.phase_duration]

                cell.dict["phenotype"].current_phase.transition_to_next_phase_args = args
                # print("_", len(cell.dict["phenotype"].current_phase.transition_to_next_phase_args))

            changed_phase, died, divides = cell.dict["phenotype"].time_step_phenotype()

            if cell.targetVolume < cell.dict["phenotype"].current_phase.volume.total:
                cell.targetVolume = cell.dict["phenotype"].current_phase.volume.total

            if changed_phase:
                cell.dict["phase_index_plus_1"] = cell.dict["phenotype"].current_phase.index + 1
                if len(self.cell_list) < 10:
                    print("@@@\nPHASE CHANGE\n@@@")

            if divides:
                if len(self.cell_list) < 10:
                    print(f"@@@\nCELL DIVISION\n@@@\ncell volume={cell.volume}")
                cells_to_divide.append(cell)

        if self.save or self.plot:
            volume_median = median(volumes)
            volume_min = min(volumes)
            volume_max = max(volumes)
            if len(time_spent_in_0):
                in_0_median = median(time_spent_in_0)
                in_0_min = min(time_spent_in_0)
                in_0_max = max(time_spent_in_0)
            else:
                in_0_median = nan
                in_0_min = nan
                in_0_max = nan
            if len(time_spent_in_1):
                in_1_median = median(time_spent_in_1)
                in_1_min = min(time_spent_in_1)
                in_1_max = max(time_spent_in_1)
            else:
                in_1_median = nan
                in_1_min = nan
                in_1_max = nan

            if self.save:
                volume_10th = quantile(volumes, 0.1)
                volume_90th = quantile(volumes, 0.9)
                volume_25th = quantile(volumes, 0.25)
                volume_75th = quantile(volumes, 0.75)

                if len(time_spent_in_0):
                    in_0_10th = quantile(time_spent_in_0, 0.1)
                    in_0_90th = quantile(time_spent_in_0, 0.9)
                    in_0_25th = quantile(time_spent_in_0, 0.25)
                    in_0_75th = quantile(time_spent_in_0, 0.75)
                else:
                    in_0_10th = nan
                    in_0_90th = nan
                    in_0_25th = nan
                    in_0_75th = nan

                if len(time_spent_in_1):
                    in_1_10th = quantile(time_spent_in_1, 0.1)
                    in_1_90th = quantile(time_spent_in_1, 0.9)
                    in_1_25th = quantile(time_spent_in_1, 0.25)
                    in_1_75th = quantile(time_spent_in_1, 0.75)
                else:
                    in_1_10th = nan
                    in_1_90th = nan
                    in_1_25th = nan
                    in_1_75th = nan

                self.volume_file.write(f"{mcs}, {volume_median}, {volume_min}, {volume_max}, "
                                       f"{volume_10th}, {volume_90th}, {volume_25th}, {volume_75th}\n")

                self.time_minus_file.write(f"{mcs}, {in_0_median}, {in_0_min}, {in_0_max}, "
                                           f"{in_0_10th}, {in_0_90th}, {in_0_25th}, {in_0_75th}\n")
                self.time_plus_file.write(f"{mcs}, {in_1_median}, {in_1_min}, {in_1_max}, "
                                          f"{in_1_10th}, {in_1_90th}, {in_1_25th}, {in_1_75th}\n")

                self.number_cells_file.write(f"{mcs}, {len(self.cell_list)}, {n_one}, {n_zero}\n")

            if self.plot:
                self.plot_win_phase.erase_all_data()
                # arguments are (name of the data series, x, y)
                self.plot_win_phase.add_data_point("N", 0, n_zero)
                self.plot_win_phase.add_data_point("N", 1, n_one)

                if not mcs % 50:
                    self.plot_win_vol.add_data_point("Median Vol", mcs, volume_median)
                    self.plot_win_vol.add_data_point("Maximum Vol", mcs, volume_max)
                    self.plot_win_vol.add_data_point("Minimum Vol", mcs, volume_min)
                    if len(time_spent_in_0):
                        self.plot_win_time.add_data_point("Median Time in Ki67-", mcs, in_0_median)

                    if len(time_spent_in_1):
                        self.plot_win_time.add_data_point("Median Time in Ki67+", mcs, median(time_spent_in_1))
                        self.plot_win_time.add_data_point("Maximum Time in Ki67+", mcs, max(time_spent_in_1))

        for cell in cells_to_divide:
            # self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume = 100  # todo: use parameter

        self.clone_parent_2_child()
        self.parent_cell.dict["phenotype"].current_phase.volume.target_cytoplasm = self.parent_cell.targetVolume
        self.parent_cell.dict["phenotype"].current_phase.volume.cytoplasm_fluid = self.parent_cell.volume
        self.parent_cell.dict["phase_index_plus_1"] = self.parent_cell.dict["phenotype"].current_phase.index + 1

        self.child_cell.dict["phenotype"].current_phase.volume.target_cytoplasm = self.parent_cell.targetVolume
        self.child_cell.dict["phenotype"].current_phase.volume.cytoplasm_fluid = self.parent_cell.targetVolume
        self.child_cell.dict["phase_index_plus_1"] = self.child_cell.dict["phenotype"].current_phase.index + 1
        if len(self.cell_list) < 10:
            print("@@@\nCHILD ATTRIBS\n@@@\n", self.child_cell.volume, self.child_cell.dict["phenotype"].time_in_phenotype,
                  self.child_cell.dict["phenotype"].current_phase,
                  self.child_cell.dict["phenotype"].current_phase.time_in_phase)
        # self.child_cell.dict["phenotype"].time_in_phenotype = 0

    def on_stop(self):
        self.finish()

    def finish(self):

        if self.save:
            self.volume_file.close()
            self.time_minus_file.close()
            self.time_plus_file.close()
            self.number_cells_file.close()