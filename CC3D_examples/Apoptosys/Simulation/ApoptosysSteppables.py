from cc3d.core.PySteppables import *

import sys

from numpy import random as rng

# sys.path.extend([abspath("../../..")])  # todo: make this more refined

# sys.path.extend(['D:\\modeling\\PhenoCellPy', 'D:/modeling/PhenoCellPy'])
sys.path.extend(['C:\\github\\PhenoCellPy', 'C:/github/PhenoCellPy'])

import Phenotypes as pheno


class ApoptosysSteppable(SteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """

        self.side = 8  # cell side
        self.dt = .5  # min/MCS

        self.to_apoptose = 5

        for cell in self.cell_list:
            cell.targetVolume = self.side*self.side
            cell.lambdaVolume = 8

    def step(self, mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """

        if mcs == 50:  # we select cells to undergo apoptosis
            self.selected_cell_ids = rng.randint(0, len(self.cell_list) + 1, self.to_apoptose, dtype=int)
            for id in self.selected_cell_ids:
                cell = self.fetch_cell_by_id(int(id))
                cell.type = self.SELECTED
                cell.lambdaVolume = 25
                apopto = pheno.phenotypes.ApoptosisStandard(dt=self.dt, nuclear_fluid=[0], nuclear_solid=[0],
                                                            cytoplasm_fluid=[.75*cell.volume],
                                                            cytoplasm_solid=[(1-.75)*cell.volume],
                                                            target_cytoplasm_to_nuclear_ratio=[0],
                                                            simulated_cell_volume=[cell.volume])

                pheno.utils.add_phenotype_to_CC3D_cell(cell, apopto)
                changed_phase, died, divides = cell.dict["phenotype"].time_step_phenotype()

        if mcs > 50:
            for id in self.selected_cell_ids:
                cell = self.fetch_cell_by_id(int(id))
                if cell is not None:
                    changed_phase, died, divides = cell.dict["phenotype"].time_step_phenotype()
                    print(cell.dict["phenotype"].current_phase.volume.cytoplasm_solid_target,
                          cell.dict["phenotype"].current_phase.volume.nuclear_solid_target,
                          cell.dict["phenotype"].current_phase.volume.target_fluid_fraction,
                          cell.dict["phenotype"].current_phase.volume.total,
                          cell.volume)
                    cell.targetVolume = cell.dict["phenotype"].current_phase.volume.total
                    cell.dict["phenotype"].current_phase.simulated_cell_volume = cell.volume

    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return
