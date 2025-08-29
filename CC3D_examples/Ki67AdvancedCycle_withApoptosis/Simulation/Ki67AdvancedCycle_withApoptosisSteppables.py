"""
BSD 3-Clause License

Copyright (c) 2023-2025, Juliano Ferrari Gianlupi and the Biocomplexity Institute, Indiana University, USA
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

from cc3d.core.PySteppables import *

from numpy import median, quantile, nan, random as rng
from numpy.random import random_sample
import math
from os.path import abspath, dirname, join

import PhenoCellPy as pcp

NUM_CELLS = 1000  # Init number of cells
CELL_WIDTH = 5  # side, pixels
dt = 1  # min/mcs
LAMBDA_VOLUME = 2.0
# APOPTOSIS_RATE  --> r1A = r2A = rqA = 1.05x10^-3 1/hr, Ghaffaraizadeh 2018 suppl 3.3.2
APOPTOSIS_RATE = 0.0000175  # 1/min --> 1.05x10^-3 1/hr * 1 hr/ 60 min

#  Q (Ki67- Quiescent phase), K1 (Ki67+ pre, prior to cell division), K2 (Ki67+ post), A (Apoptosis)

# check whether to transition from Ki67+ pre to mitosis and Ki67+ post:
def Ki67pos_transition(*args) -> bool:
    # print(len(args), print(args))
    # args = [cc3d cell volume, phase's target volume, time in phase, phase duration
    # print("Ki67pos_transition(): ", args[0], ", ", args[1],",time in phase: ", args[2], ", Phase duration: ", args[3])
    return args[0] >= args[1] and args[2] > args[3]


class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.track_cell_level_scalar_attribute(field_name='phase_index_plus_1', attribute_name='phase_index_plus_1')
        self.target_volume = None
        self.doubling_volume = None
        self.volume_conversion_unit = None


    def start(self):
        side = CELL_WIDTH
        self.target_volume = side * side
        self.doubling_volume = 2 * self.target_volume

        cell_area = math.pi * side * side
        total_blob_area = cell_area * NUM_CELLS
        current_radius = math.sqrt(total_blob_area / math.pi)
        if current_radius < self.dim.x / 2:
            blob_radius = int(current_radius * 0.57) # fudge factor to produce a full circle of cells
        else:
            blob_radius = self.dim.x // 35  # in lattice units

        center_x = self.dim.x // 2 - side // 2
        center_y = self.dim.y // 2 - side // 2

        # Define the bounding box around the blob
        min_x = max(0, center_x - blob_radius)
        max_x = min(self.dim.x, center_x + blob_radius)
        min_y = max(0, center_y - blob_radius)
        max_y = min(self.dim.y, center_y + blob_radius)

        created_cells = set()

        for x in range(min_x, max_x, side):
            for y in range(min_y, max_y, side):
                dx = x - center_x
                dy = y - center_y
                distance = math.sqrt(dx**2 + dy**2)
                if distance <= blob_radius:
                    if (x, y) not in created_cells and len(created_cells) <= NUM_CELLS:
                        cell = self.new_cell(self.CELL)
                        cell.type = self.CELL
                        self.cell_field[x:x + side - 1, y: y + side -1, 0] = cell
                        created_cells.add((x, y))

        ki67_advanced_modified_transition = pcp.phenotypes.Ki67Advanced(dt=dt,
                                                        nuclear_volume_change_rate=[None, 0.0055, None],
                                                        cytoplasm_volume_change_rate=[None, 0.0045, None],
                                                        fluid_change_rate=[None, 0.05, None],
                                                        phase_durations=[74.35 * 60, 13 * 60, 2.5 * 60],
                                                        check_transition_to_next_phase_functions=
                                                        [None, Ki67pos_transition, None],
                                                        check_transition_to_next_phase_functions_args=[None,
                                                                                        [-9, 1, -9, 1], None])

        self.volume_conversion_unit = self.target_volume / ki67_advanced_modified_transition.current_phase.volume.total

        for cell in self.cell_list:
            cell.targetVolume = self.target_volume
            cell.lambdaVolume = LAMBDA_VOLUME # 2.0
            pcp.utils.add_phenotype_to_CC3D_cell(cell, ki67_advanced_modified_transition)
            cell.dict["phase_index_plus_1"] = cell.dict["phenotype"].current_phase.index + 1

        self.shared_steppable_vars["constraints"] = self


class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)

        self.constraint_vars = None
        self.previous_number_cells = 0

        self.plot = True
        self.save = False

        if self.save:
            self.save_loc = dirname(abspath(__file__))

            self.volume_file = open(join(self.save_loc, "volume.dat"), "w+")
            self.volume_file.write("MCS, Median, Min, Max, 10th, 90th, 25th, 75th\n")

            self.time_minus_file = open(join(self.save_loc, "time_in_Ki67-.dat"), "w+")
            self.time_minus_file.write("MCS, Median, Min, Max, 10th, 90th, 25th, 75th\n")

            self.time_plus_file = open(join(self.save_loc, "time_in_Ki67+pre.dat"), "w+")
            self.time_plus_file.write("MCS, Median, Min, Max, 10th, 90th, 25th, 75th\n")

            self.time_plus_post_file = open(join(self.save_loc, "time_in_Ki67+post.dat"), "w+")
            self.time_plus_post_file.write("MCS, Median, Min, Max, 10th, 90th, 25th, 75th\n")

            self.ratio_cells_in_phase_file = open(join(self.save_loc, "ratio_cells_in_phase.dat"), "w+")
            self.ratio_cells_in_phase_file.write("hr, Ki67-, Ki67+pre, Ki67+post, Apoptosis\n")

            self.number_cells_file = open(join(self.save_loc, "number_cells.dat"), "w+")
            self.number_cells_file.write("MCS, N tot, N Apop, N+ post, N+ pre, N-\n")

    def start(self):
        self.constraint_vars = self.shared_steppable_vars["constraints"]
        if self.plot:
            self.plot_win_vol = self.add_new_plot_window(title='Volume metrics',
                                                         x_axis_title='MonteCarlo Step (MCS)',
                                                         y_axis_title='Volume', x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=True,
                                                         config_options={"legend": True})

            self.plot_win_vol.add_plot("Median Vol", style='Lines', color='yellow', size=5)
            self.plot_win_vol.add_plot("Minimum Vol", style='Lines', color='blue', size=5)
            self.plot_win_vol.add_plot("Maximum Vol", style='Lines', color='red', size=5)

            self.plot_win_time = self.add_new_plot_window(title='Time spent in each phase',
                                                          x_axis_title='MonteCarlo Step (MCS)',
                                                          y_axis_title='Total MCSs', x_scale_type='linear',
                                                          y_scale_type='linear',
                                                          grid=True,
                                                          config_options={"legend": True})

            #self.plot_win_time.add_plot("Maximum Time in Ki67-", style='Lines', color='blue', size=5)
            self.plot_win_time.add_plot("Median Time in Ki67-", style='Lines', color='yellow', size=5)
            self.plot_win_time.add_plot("Median Time in Ki67+ pre", style='Lines', color='red', size=5)
            self.plot_win_time.add_plot("Median Time in Ki67+ post", style='Lines', color='white', size=5)
            self.plot_win_time.add_plot("Median Time in Apoptosis", style='Lines', color='green', size=5)

            self.plot_win_phase_ratio = self.add_new_plot_window(title='Ratio of cells in each phase.',
                                                    x_axis_title='hours (MCS * dt(min/mcs) / 60)',
                                                    y_axis_title='# cells in phase / total cells', x_scale_type='linear',
                                                    y_scale_type='linear',
                                                    grid=True,
                                                    config_options={"legend": True})

            self.plot_win_phase_ratio.add_plot("Ratio of cells in phase Ki67-", style='Lines', color='yellow', size=5)
            self.plot_win_phase_ratio.add_plot("Ratio of cells in phase Ki67+ pre", style='Lines', color='red', size=5)
            self.plot_win_phase_ratio.add_plot("Ratio of cells in phase Ki67+ post", style='Lines', color='white', size=5)
            self.plot_win_phase_ratio.add_plot("Ratio of cells in phase Apoptosis", style='Lines', color='green', size=5)

            self.plot_win_number = self.add_new_plot_window(title='Number of cells',
                                                            x_axis_title='MonteCarlo Step (MCS)',
                                                            y_axis_title='# of cells', x_scale_type='linear',
                                                            y_scale_type='linear',
                                                            grid=True)
            self.plot_win_number.add_plot("N", style='Lines', color='red', size=5)

            self.plot_win_phase = self.add_new_plot_window(title='Number of cells in each phase',
                                                           x_axis_title='Phase: 0, 1, 2, 3 (Apop)',
                                                           y_axis_title='# cells', x_scale_type='linear',
                                                           y_scale_type='linear',
                                                           grid=True)

            self.plot_win_phase.add_plot("N", style='Points', color='red', size=20)

    def step(self, mcs):

        if not mcs and self.plot:
            self.plot_win_number.add_data_point("N", mcs, len(self.cell_list))
            self.previous_number_cells = len(self.cell_list)

        elif not mcs % 50 and self.plot:
            self.plot_win_number.add_data_point("N", mcs, len(self.cell_list))
            self.previous_number_cells = len(self.cell_list)

        cells_to_divide = []

        n_zero = 0
        n_one = 0
        n_two = 0
        n_apop = 0

        volumes = []
        # Time spent in Phase 0 (TQ): 4.59 * 60 min, Phase 1 (T1): 13.0 * 60.0, Phase 2 (T2): 2.5 * 60  (phenotypes.py: Ki67Advanced)
        #               Phase 3 (TA) Apoptosis: 8.6 * 60
        # From Ghaffarizadeh 2018 suppl sec 4.:
        # For both main examples, we used the Ki-67 advanced model (legacy variant), with T1 = 13 hours, T2 = 2.5 hours,
        # ⟨T⟩Q = 74.35 hours, and TA = 8.6 hours.

        time_spent_in_0 = []
        time_spent_in_1 = []
        time_spent_in_2 = []
        time_spent_in_apop = []
        total_cells = len(self.cell_list)

        for cell in self.cell_list:
            if cell.dict["phenotype"].name != "Standard apoptosis model":
                volumes.append(cell.volume)
                cell.dict["phenotype"].current_phase.simulated_cell_volume = cell.volume
                if cell.dict["phenotype"].current_phase.index == 0:
                    n_zero += 1
                    time_spent_in_0.append(cell.dict["phenotype"].current_phase.time_in_phase)
                    #if cell.dict["phenotype"].current_phase.time_in_phase > expected_time_in_0 / dt:
                        #print("Phase 0 : in longer than max: ", cell.dict["phenotype"].current_phase.time_in_phase)
                elif cell.dict["phenotype"].current_phase.index == 1:
                    n_one += 1
                    time_spent_in_1.append(cell.dict["phenotype"].current_phase.time_in_phase)
                    # args = [cc3d cell volume, doubling volume, time in phase, phase duration
                    args = [
                        cell.volume,
                        .9 * self.constraint_vars.doubling_volume,  # we use 90% of the doubling volume because cc3d cells
                        # will always be slightly below their target due to the contact energy
                        cell.dict["phenotype"].current_phase.time_in_phase + cell.dict["phenotype"].dt,
                        cell.dict["phenotype"].current_phase.phase_duration]

                    cell.dict["phenotype"].current_phase.check_transition_to_next_phase_function_args = args
                elif cell.dict["phenotype"].current_phase.index == 2:
                    n_two += 1
                    time_spent_in_2.append(cell.dict["phenotype"].current_phase.time_in_phase)

                changed_phase, should_be_removed, divides = cell.dict["phenotype"].time_step_phenotype()
                converted_volume = self.constraint_vars.volume_conversion_unit * \
                                    cell.dict["phenotype"].current_phase.volume.total
                cell.targetVolume = converted_volume

                if changed_phase:
                    cell.dict["phase_index_plus_1"] = cell.dict["phenotype"].current_phase.index + 1
                    print("@@@\nPHASE CHANGE\n@@@")

                if divides:
                    print(f"@@@\nCELL DIVISION\n@@@\ncell volume={cell.volume}")
                    cells_to_divide.append(cell)
            else:  # Apoptosis phenotype:
                n_apop += 1
                time_spent_in_apop.append(cell.dict["phenotype"].current_phase.time_in_phase)
                print(f"@@@\nCELL DEATH\n@@@\ncell type={cell.type}")

        if total_cells < 1:
            total_cells = 1
        cell_in_phase_0_ratio = (n_zero / total_cells)
        cell_in_phase_1_ratio = (n_one / total_cells)
        cell_in_phase_2_ratio = (n_two / total_cells)
        cell_in_phase_ap_ratio = (n_apop / total_cells)

        if self.save or self.plot:
            time_hour = dt * mcs / 60
            if len(volumes) > 0:
                volume_median = median(volumes)
                volume_min = min(volumes)
                volume_max = max(volumes)
            else:
                volume_median = nan
                volume_min = nan
                volume_max = nan
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
            if len(time_spent_in_2):
                in_2_median = median(time_spent_in_2)
                in_2_min = min(time_spent_in_2)
                in_2_max = max(time_spent_in_2)
            else:
                in_2_median = nan
                in_2_min = nan
                in_2_max = nan
            if len(time_spent_in_apop):
                in_apop_median = median(time_spent_in_apop)
                in_apop_min = min(time_spent_in_apop)
                in_apop_max = max(time_spent_in_apop)
            else:
                in_apop_median = nan
                in_apop_min = nan
                in_apop_max = nan

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

                if len(time_spent_in_2):
                    in_2_10th = quantile(time_spent_in_2, 0.1)
                    in_2_90th = quantile(time_spent_in_2, 0.9)
                    in_2_25th = quantile(time_spent_in_2, 0.25)
                    in_2_75th = quantile(time_spent_in_2, 0.75)
                else:
                    in_2_10th = nan
                    in_2_90th = nan
                    in_2_25th = nan
                    in_2_75th = nan

                if len(time_spent_in_apop):
                    in_apop_10th = quantile(time_spent_in_apop, 0.1)
                    in_apop_90th = quantile(time_spent_in_apop, 0.9)
                    in_apop_25th = quantile(time_spent_in_apop, 0.25)
                    in_apop_75th = quantile(time_spent_in_apop, 0.75)
                else:
                    in_apop_10th = nan
                    in_apop_90th = nan
                    in_apop_25th = nan
                    in_apop_75th = nan

                self.volume_file.write(f"{mcs}, {volume_median}, {volume_min}, {volume_max}, "
                                       f"{volume_10th}, {volume_90th}, {volume_25th}, {volume_75th}\n")

                self.time_minus_file.write(f"{mcs}, {in_0_median}, {in_0_min}, {in_0_max}, "
                                           f"{in_0_10th}, {in_0_90th}, {in_0_25th}, {in_0_75th}\n")
                self.time_plus_file.write(f"{mcs}, {in_1_median}, {in_1_min}, {in_1_max}, "
                                          f"{in_1_10th}, {in_1_90th}, {in_1_25th}, {in_1_75th}\n")
                self.time_plus_post_file.write(f"{mcs}, {in_2_median}, {in_2_min}, {in_2_max}, "
                                               f"{in_2_10th}, {in_2_90th}, {in_2_25th}, {in_2_75th}\n")
                self.time_plus_file.write(f"{mcs}, {in_apop_median}, {in_apop_min}, {in_apop_max}, "
                                          f"{in_apop_10th}, {in_apop_90th}, {in_apop_25th}, {in_apop_75th}\n")

                self.ratio_cells_in_phase_file.write(f"{time_hour}, {cell_in_phase_0_ratio}, {cell_in_phase_1_ratio}, "
                                                     f"{cell_in_phase_2_ratio}, {cell_in_phase_ap_ratio}\n")
                self.number_cells_file.write(f"{mcs}, {len(self.cell_list)}, {n_apop}, {n_two}, {n_one}, {n_zero}\n")

            if self.plot:
                self.plot_win_phase.erase_all_data()
                # arguments are (name of the data series, x, y)
                self.plot_win_phase.add_data_point("N", 0, n_zero)
                self.plot_win_phase.add_data_point("N", 1, n_one)
                self.plot_win_phase.add_data_point("N", 2, n_two)
                self.plot_win_phase.add_data_point("N", 3, n_apop)

                if not mcs % 50:
                    self.plot_win_vol.add_data_point("Median Vol", mcs, volume_median)
                    self.plot_win_vol.add_data_point("Maximum Vol", mcs, volume_max)
                    self.plot_win_vol.add_data_point("Minimum Vol", mcs, volume_min)
                    if len(time_spent_in_0):
                        #self.plot_win_time.add_data_point("Maximum Time in Ki67-", mcs, max(time_spent_in_0))
                        self.plot_win_time.add_data_point("Median Time in Ki67-", mcs, median(time_spent_in_0))

                    if len(time_spent_in_1):
                        self.plot_win_time.add_data_point("Median Time in Ki67+ pre", mcs, in_1_median)

                    if len(time_spent_in_2):
                        self.plot_win_time.add_data_point("Median Time in Ki67+ post", mcs, in_2_median)
                    if len(time_spent_in_apop):
                        self.plot_win_time.add_data_point("Median Time in Apoptosis", mcs, in_apop_median)

                    self.plot_win_phase_ratio.add_data_point("Ratio of cells in phase Ki67-", time_hour, cell_in_phase_0_ratio)
                    self.plot_win_phase_ratio.add_data_point("Ratio of cells in phase Ki67+ pre", time_hour, cell_in_phase_1_ratio)
                    self.plot_win_phase_ratio.add_data_point("Ratio of cells in phase Ki67+ post", time_hour, cell_in_phase_2_ratio)
                    self.plot_win_phase_ratio.add_data_point("Ratio of cells in phase Apoptosis", time_hour, cell_in_phase_ap_ratio)

        for cell in cells_to_divide:
            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # resetting target volume
        self.parent_cell.dict["phenotype"].current_phase.volume.nuclear_solid = self.parent_cell.dict["phenotype"].current_phase.volume.nuclear_solid / 2
        self.parent_cell.dict["phenotype"].current_phase.volume.cytoplasm_solid = self.parent_cell.dict["phenotype"].current_phase.volume.cytoplasm_solid / 2
        converted_volume = self.constraint_vars.volume_conversion_unit * \
                           self.parent_cell.dict["phenotype"].current_phase.volume.total
        self.parent_cell.targetVolume = converted_volume

        self.clone_parent_2_child()
        self.child_cell.dict["phenotype"] = self.parent_cell.dict["phenotype"].copy()
        self.parent_cell.dict["phase_index_plus_1"] = self.parent_cell.dict["phenotype"].current_phase.index + 1

        self.child_cell.dict["phase_index_plus_1"] = self.child_cell.dict["phenotype"].current_phase.index + 1
        self.child_cell.dict["phenotype"].time_in_phenotype = 0
        if len(self.cell_list) < CELL_FIELD * CELL_FIELD:
            print("@@@\nCHILD ATTRIBS\n@@@\n", self.child_cell.volume,
                  self.child_cell.dict["phenotype"].time_in_phenotype,
                  self.child_cell.dict["phenotype"].current_phase,
                  self.child_cell.dict["phenotype"].current_phase.time_in_phase)

    def on_stop(self):
        self.finish()

    def finish(self):

        if self.save:
            self.volume_file.close()
            self.time_minus_file.close()
            self.time_plus_file.close()
            self.ratio_cells_in_phase_file.close()
            self.number_cells_file.close()


class ApoptosisSteppable(SteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.rate_apop = APOPTOSIS_RATE  # 1/min
        #print("apoptosis: ", float(self.rate_apop * dt))

        self.apopto = pcp.phenotypes.ApoptosisStandard(dt=dt, )
        self.selected_cell_ids = []  # list of cells being destroyed

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        self.target_volume = CELL_WIDTH * CELL_WIDTH
        self.volume_conversion_unit = self.target_volume / self.apopto.current_phase.volume.total

    def step(self, mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        if random_sample() < float(self.rate_apop * dt):
            if len(self.cell_list) > 0:
                cid = rng.randint(1, len(self.cell_list) + 1)  # no cell id of '0' ?
                if type(self.fetch_cell_by_id(int(cid))) == "None":
                    print(" cell with id ", cid, " has no type associated with it!")
                try:
                    if self.fetch_cell_by_id(int(cid)).type != self.APOPTOSIS:
                        cell = self.fetch_cell_by_id(int(cid))
                        cell.type = self.APOPTOSIS
                        cell.dict["phase_index_plus_1"] = 4  # add to cell tracking scaler attribute defined earlier.
                                                             # Assumes only 3 other phases for cell.
                        cell.lambdaVolume = 25  # the cell is not alive anymore, so it should be pretty stiff as it can't
                        # reshape itself actively.
                        cell.dict["phenotype"] = None  # Remove any existing phenotype
                        self.selected_cell_ids.append(int(cid))
                        pcp.utils.add_phenotype_to_CC3D_cell(cell, self.apopto)
                        changed_phase, should_be_removed, divides = cell.dict["phenotype"].time_step_phenotype()
                        if should_be_removed:
                            self.delete_cell(cell)
                except AttributeError:
                    print(" Cell id: ", cid, ", is of type 'NoneType' or has no cell type defined")

            # reduce volumes of cells until they are removed
            for cid in self.selected_cell_ids:
                cell = self.fetch_cell_by_id(int(cid))
                if cell is not None:
                    changed_phase, should_be_removed, divides = cell.dict["phenotype"].time_step_phenotype()
                    # print(cell.dict["phenotype"].current_phase.volume.cytoplasm_solid_target,
                          #  cell.dict["phenotype"].current_phase.volume.nuclear_solid_target,
                          #  cell.dict["phenotype"].current_phase.volume.target_fluid_fraction,
                          #  cell.dict["phenotype"].current_phase.volume.total,
                          #  cell.volume)
                    cell.targetVolume = self.volume_conversion_unit * \
                                        cell.dict["phenotype"].current_phase.volume.total
                    cell.dict["phenotype"].current_phase.simulated_cell_volume = cell.volume
                    if should_be_removed:
                        print(" ----- > cell type removed (type 2 is correct one):", cell.type)
                        self.delete_cell(cell)

    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return
