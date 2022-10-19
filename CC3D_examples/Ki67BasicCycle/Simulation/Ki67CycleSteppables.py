from cc3d.cpp.PlayerPython import *
from cc3d import CompuCellSetup

from cc3d.core.PySteppables import *

from numpy import median, quantile, nan

import sys
from os.path import abspath, dirname, join

# sys.path.extend([abspath("../../..")])  # todo: make this more refined

sys.path.extend(['D:\\modeling\\PhenoCellPy', 'D:/modeling/PhenoCellPy'])

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

        ki67_basic = pheno.phenotypes.Ki67Basic(dt=dt, target_volumes=[side * side, side * side],
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

        self.save = True

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

        time_spent_in_0 = []
        time_spent_in_1 = []

        for cell in self.cell_list:
            volumes.append(cell.volume)
            if cell.dict["cycle"].current_phase.index == 0:
                n_zero += 1
                time_spent_in_0.append(cell.dict["cycle"].current_phase.time_in_phase)
            elif cell.dict["cycle"].current_phase.index == 1:
                n_one += 1
                time_spent_in_1.append(cell.dict["cycle"].current_phase.time_in_phase)
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

    def on_stop(self):
        self.finish()

    def finish(self):

        if self.save:
            self.volume_file.close()
            self.time_minus_file.close()
            self.time_plus_file.close()
            self.number_cells_file.close()
