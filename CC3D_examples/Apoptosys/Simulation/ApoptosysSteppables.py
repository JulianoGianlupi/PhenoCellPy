"""
BSD 3-Clause License

Copyright (c) 2023, Juliano Ferrari Gianlupi
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

import sys

from numpy import random as rng

# sys.path.extend([abspath("../../..")])  # todo: make this more refined

# sys.path.extend(['D:\\modeling\\PhenoCellPy', 'D:/modeling/PhenoCellPy'])
sys.path.extend(['C:\\github\\PhenoCellPy', 'C:/github/PhenoCellPy'])

import PhenoCellPy as pheno


class ApoptosysSteppable(SteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """

        self.side = 8  # cell side
        self.dt = .1  # min/MCS

        self.to_apoptose = 5

        self.target_volume = self.side * self.side

        self.apopto = pheno.phenotypes.ApoptosisStandard(dt=self.dt)

        self.volume_conversion_unit = self.target_volume / self.apopto.current_phase.volume.total

        for cell in self.cell_list:
            cell.targetVolume = self.target_volume
            cell.lambdaVolume = 8

    def step(self, mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """

        if mcs == 50:  # we select cells to undergo apoptosis
            self.selected_cell_ids = rng.randint(0, len(self.cell_list) + 1, self.to_apoptose)
            for cid in self.selected_cell_ids:
                cell = self.fetch_cell_by_id(int(cid))
                cell.type = self.SELECTED
                cell.lambdaVolume = 25  # the cell is not alive anymore, so it should be pretty stiff as it can't
                # reshape itself actively

                pheno.utils.add_phenotype_to_CC3D_cell(cell, self.apopto)
                changed_phase, should_be_removed, divides = cell.dict["phenotype"].time_step_phenotype()
                if should_be_removed:
                    self.delete_cell(cell)

        if mcs > 50:
            for cid in self.selected_cell_ids:
                cell = self.fetch_cell_by_id(int(cid))
                if cell is not None:
                    changed_phase, should_be_removed, divides = cell.dict["phenotype"].time_step_phenotype()
                    print(cell.dict["phenotype"].current_phase.volume.cytoplasm_solid_target,
                          cell.dict["phenotype"].current_phase.volume.nuclear_solid_target,
                          cell.dict["phenotype"].current_phase.volume.target_fluid_fraction,
                          cell.dict["phenotype"].current_phase.volume.total,
                          cell.volume)
                    cell.targetVolume = self.volume_conversion_unit * \
                                        cell.dict["phenotype"].current_phase.volume.total
                    cell.dict["phenotype"].current_phase.simulated_cell_volume = cell.volume
                    if should_be_removed:
                        self.delete_cell(cell)

    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return
