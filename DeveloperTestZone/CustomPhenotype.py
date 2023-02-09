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

import tissue_forge as tf
import numpy as np

import sys
from os.path import abspath
import time

sys.path.extend([abspath("../")])  # todo: make this more refined

import Phenotypes as pheno


def get_radius_sphere(volume):
    return ((1 / (np.pi * 4 / 3)) * volume) ** (1 / 3)


# potential cutoff distance
cutoff = 3

# space set up
dim = [50, 50, 50]
tf.init(dim=dim, cutoff=cutoff)

pot = tf.Potential.morse(d=3, a=5, min=-0.8, max=2)

# Particle types

# volume = 10

mass = 40

radius = .4

global density
density = mass / ((4 / 3) * np.pi * radius * radius * radius)

dt = 10  # min/time step

stable_phase_0 = pheno.phases.Phase(index=0, previous_phase_index=-1, next_phase_index=1, dt=dt,
                                    time_unit="min", space_unit="micrometer", name="stable0",
                                    division_at_phase_exit=False, removal_at_phase_exit=False, fixed_duration=True,
                                    phase_duration=20, entry_function=None, exit_function=None, arrest_function=None,
                                    check_transition_to_next_phase_function=None,
                                    simulated_cell_volume=np.pi * (4 / 3) * radius * radius * radius,
                                    cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None,
                                    calcification_rate=None, target_fluid_fraction=None, nuclear_fluid=None,
                                    nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                                    cytoplasm_solid=None, cytoplasm_solid_target=None,
                                    target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None,
                                    fluid_change_rate=None, relative_rupture_volume=None,
                                    user_phase_time_step=None, user_phase_time_step_args=(None,))


def grow_phase_transition(*args):
    return args[0] >= args[1] and args[2] > args[4]

def double_target_volumes(self, *none):
    self.volume.nuclear_solid_target *= 2
    self.volume.cytoplasm_solid_target *= 2


grow_phase = pheno.phases.Phase(index=1, previous_phase_index=0, next_phase_index=2, dt=dt,
                                time_unit="min", space_unit="micrometer", name="grow",
                                division_at_phase_exit=False, removal_at_phase_exit=False, fixed_duration=True,
                                phase_duration=120, entry_function=double_target_volumes,
                                entry_function_args=[None],
                                exit_function=None, arrest_function=None,
                                check_transition_to_next_phase_function=grow_phase_transition,
                                check_transition_to_next_phase_function_args=[0, 9, 0, 9],
                                simulated_cell_volume=np.pi * (4 / 3) * radius * radius * radius,
                                cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None,
                                calcification_rate=None, target_fluid_fraction=None, nuclear_fluid=None,
                                nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                                cytoplasm_solid=None, cytoplasm_solid_target=None,
                                target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None,
                                fluid_change_rate=None, relative_rupture_volume=None,
                                user_phase_time_step=None, user_phase_time_step_args=(None,))

stable_phase_1 = pheno.phases.Phase(index=2, previous_phase_index=1, next_phase_index=3, dt=dt,
                                    time_unit="min", space_unit="micrometer", name="stable1",
                                    division_at_phase_exit=False, removal_at_phase_exit=False, fixed_duration=True,
                                    phase_duration=30, entry_function=None, exit_function=None, arrest_function=None,
                                    check_transition_to_next_phase_function=None,
                                    simulated_cell_volume=np.pi * (4 / 3) * radius * radius * radius,
                                    cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None,
                                    calcification_rate=None, target_fluid_fraction=None, nuclear_fluid=None,
                                    nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                                    cytoplasm_solid=None, cytoplasm_solid_target=None,
                                    target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None,
                                    fluid_change_rate=None, relative_rupture_volume=None,
                                    user_phase_time_step=None, user_phase_time_step_args=(None,))


def shrink_phase_transition(self, *none):
    time_check = np.random.uniform() < (1 - np.exp(-self.dt / self.phase_duration))
    volume_check = self.volume.total <= 1.1 * self.volume.total_target
    return time_check and volume_check


shrink_phase = pheno.phases.Phase(index=3, previous_phase_index=2, next_phase_index=0, dt=dt,
                                  time_unit="min", space_unit="micrometer", name="shrink",
                                  division_at_phase_exit=False, removal_at_phase_exit=False, fixed_duration=False,
                                  phase_duration=60, entry_function=pheno.phases.Phase._halve_target_volume,
                                  entry_function_args=[None], exit_function=None, arrest_function=None,
                                  check_transition_to_next_phase_function=shrink_phase_transition,
                                  check_transition_to_next_phase_function_args=[None],
                                  simulated_cell_volume=np.pi * (4 / 3) * radius * radius * radius,
                                  cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None,
                                  calcification_rate=None, target_fluid_fraction=None, nuclear_fluid=None,
                                  nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                                  cytoplasm_solid=None, cytoplasm_solid_target=None,
                                  target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None,
                                  fluid_change_rate=None, relative_rupture_volume=None,
                                  user_phase_time_step=None, user_phase_time_step_args=(None,))


custom_pheno = pheno.phenotypes.Phenotype(name="oscillate volume with rests", dt=dt, time_unit="min",
                                          space_unit="micrometer", phases=[stable_phase_0, grow_phase, stable_phase_1,
                                                                           shrink_phase],
                                          quiescent_phase=False, starting_phase_index=0, user_phenotype_time_step=None,
                                          user_phenotype_time_step_args=[None,])

global volume_conversion_unit
volume_conversion_unit = mass / custom_pheno.current_phase.volume.total


class CellType(tf.ParticleTypeSpec):
    mass = mass
    target_temperature = 0
    radius = radius
    dynamics = tf.Overdamped
    cycle = custom_pheno


Cell = CellType.get()

tf.bind.types(pot, Cell, Cell)

rforce = tf.Force.random(mean=0, std=50)

# bind it just like any other force
tf.bind.force(rforce, Cell)

first_cell = Cell([d // 2 for d in dim])
first_cell.cycle = custom_pheno

global cells_cycles

cells_cycles = {f"{first_cell.id}": custom_pheno}


def step_cycle_and_divide(event):
    for p in Cell.items():
        pcycle = cells_cycles[f"{p.id}"]
        pcycle.current_phase.simulated_cell_volume = p.mass * density
        if pcycle.current_phase.name == "grow" or pcycle.current_phase.name == "shrink":
            phase_change, should_be_removed, division = pcycle.time_step_phenotype()
        else:
            phase_change, should_be_removed, division = pcycle.time_step_phenotype()
        print(p.id, p.radius, pcycle.current_phase.name, pcycle.current_phase.volume.total)

        if phase_change and len(Cell.items()) < 10:
            print("@@@\nPHASE CHANGE\n@@@")
            # time.sleep(1)

        radius = get_radius_sphere(volume_conversion_unit * pcycle.current_phase.volume.total)

        # book-keeping, making sure the simulated cell grows
        # if p.radius < radius:
        p.radius = radius
        p.mass = ((4 / 3) * np.pi * radius * radius * radius) * density

        # if division occurs, divide
        if division:
            print("@@@\nDIVISION\n@@@")
            # save cell attribs to halve later
            cur_mass = p.mass

            # divide and reasign attribs (is this step necessary?)
            child = p.split()

            cells_cycles[f"{child.id}"] = custom_pheno

            child.mass = p.mass = cur_mass / 2
            child.radius = p.radius = get_radius_sphere((cur_mass / 2) / density)

            cells_cycles[f"{child.id}"].volume = child.mass * density
            cells_cycles[f"{child.id}"].simulated_cell_volume = child.mass * density

    return 0


tf.event.on_time(invoke_method=step_cycle_and_divide, period=.9 * tf.Universe.dt)

# run the simulator interactive
tf.run()
