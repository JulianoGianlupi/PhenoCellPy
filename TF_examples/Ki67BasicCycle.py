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

ki67_basic = pheno.phenotypes.Ki67Basic(dt=dt)

global volume_conversion_unit
volume_conversion_unit = (4 / 3) * np.pi * radius * radius * radius/ki67_basic.current_phase.volume.total


class CellType(tf.ParticleTypeSpec):
    mass = mass
    target_temperature = 0
    radius = radius
    dynamics = tf.Overdamped
    cycle = ki67_basic


Cell = CellType.get()

tf.bind.types(pot, Cell, Cell)

rforce = tf.Force.random(mean=0, std=50)

# bind it just like any other force
tf.bind.force(rforce, Cell)

first_cell = Cell([d // 2 for d in dim])
first_cell.cycle = ki67_basic

global cells_cycles

cells_cycles = {f"{first_cell.id}": ki67_basic}


def step_cycle_and_divide(event):
    for p in Cell.items():
        pcycle = cells_cycles[f"{p.id}"]
        pcycle.current_phase.simulated_cell_volume = p.mass * density
        phase_change, should_be_removed, division = pcycle.time_step_phenotype()

        if phase_change and len(Cell.items()) < 10:
            print("@@@\nPHASE CHANGE\n@@@")
            # time.sleep(1)

        radius = get_radius_sphere(volume_conversion_unit*pcycle.current_phase.volume.total)

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

            cells_cycles[f"{child.id}"] = ki67_basic

            child.mass = p.mass = cur_mass / 2
            child.radius = p.radius = get_radius_sphere((cur_mass / 2) / density)

            cells_cycles[f"{child.id}"].volume = child.mass * density
            cells_cycles[f"{child.id}"].simulated_cell_volume = child.mass * density

    return 0


tf.event.on_time(invoke_method=step_cycle_and_divide, period=.9*tf.Universe.dt)

# run the simulator interactive
tf.run()
