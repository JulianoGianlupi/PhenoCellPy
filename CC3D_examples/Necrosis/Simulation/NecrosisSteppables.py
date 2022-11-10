from cc3d.core.PySteppables import *

import sys

from numpy import random as rng

# sys.path.extend([abspath("../../..")])  # todo: make this more refined

# sys.path.extend(['D:\\modeling\\PhenoCellPy', 'D:/modeling/PhenoCellPy'])
sys.path.extend(['C:\\github\\PhenoCellPy', 'C:/github/PhenoCellPy'])

import Phenotypes as pheno


class NecrosisSteppable(SteppableBasePy):

    def __init__(self, frequency=1):

        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """
        Called before MCS=0 while building the initial simulation
        """

        self.side = 8  # cell side
        self.dt = 2.5  # min/MCS

        self.to_necrose = 5  # how many cells we will necrose

        for cell in self.cell_list:
            cell.targetVolume = self.side * self.side
            cell.lambdaVolume = 8
            necrotic_phenotype = pheno.phenotypes.NecrosisStandard(dt=self.dt,
                                                                   nuclear_fluid=[0, 0],
                                                                   nuclear_solid=[0, 0],
                                                                   cytoplasm_fluid=[.75 * cell.volume,
                                                                                    .75 * cell.volume],
                                                                   # the default volume
                                                                   # model uses a .75 fluid fraction
                                                                   cytoplasm_solid=[(1 - .75) * cell.volume,
                                                                                    (1 - .75) * cell.volume],
                                                                   simulated_cell_volume=cell.volume)

    def step(self, mcs):
        """
        Called every frequency MCS while executing the simulation
        
        :param mcs: current Monte Carlo step
        """
        if mcs == 50:  # we select cells to undergo necrosis
            self.selected_cell_ids = rng.randint(0, len(self.cell_list) + 1, self.to_necrose)
            for cid in self.selected_cell_ids:
                cell = self.fetch_cell_by_id(int(cid))
                cell.type = self.NECROTIC
                cell.lambdaVolume = 25  # the cell is not alive anymore, so it should be pretty stiff as it can't
                # reshape itself actively

                necrotic_phenotype = pheno.phenotypes.NecrosisStandard(dt=self.dt,
                                                                       nuclear_fluid=[0, 0],
                                                                       nuclear_solid=[0, 0],
                                                                       cytoplasm_fluid=[.75 * cell.volume,
                                                                                        .75 * cell.volume],
                                                                       # the default volume
                                                                       # model uses a .75 fluid fraction
                                                                       cytoplasm_solid=[(1 - .75) * cell.volume,
                                                                                        (1 - .75) * cell.volume],
                                                                       simulated_cell_volume=cell.volume)
                pheno.utils.add_phenotype_to_CC3D_cell(cell, necrotic_phenotype)

        if mcs > 50:
            for cid in self.selected_cell_ids:
                cell = self.fetch_cell_by_id(int(cid))
                if cell is not None:  # if the cell has died (disappeared, deleted from the simulation) cell is a null
                    # object
                    # todo: monitor the phase change, change type when it happens, etc
                    changed_phase, died, divides = cell.dict["phenotype"].time_step_phenotype()
                    print("TARGET VOLUMES:\n",
                          "Cyto solid:", cell.dict["phenotype"].current_phase.volume.cytoplasm_solid_target,
                          ". Nucl solid:", cell.dict["phenotype"].current_phase.volume.nuclear_solid_target,
                          ". Fluid frac solid:", cell.dict["phenotype"].current_phase.volume.target_fluid_fraction,
                          "\n CURRENT VOLUMES:\n",
                          "TOTAL:", cell.dict["phenotype"].current_phase.volume.total,
                          "cyto:", cell.dict["phenotype"].current_phase.volume.cytoplasm,
                          ". nucl:", cell.dict["phenotype"].current_phase.volume.nuclear,
                          ". fluid:", cell.dict["phenotype"].current_phase.volume.fluid,
                          "\n SIMULATED CELL:", cell.volume)
                    cell.targetVolume = cell.dict["phenotype"].current_phase.volume.total
                    cell.dict["phenotype"].current_phase.simulated_cell_volume = cell.volume
                    if changed_phase:
                        print("CELL BURST!")
                        if cell.dict["phenotype"].current_phase.name == "Necrotic (lysed)":
                            # if the cell has ruptured
                            cell.type = self.RUPTURED

    def finish(self):
        """
        Called after the last MCS to wrap up the simulation
        """

    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """
