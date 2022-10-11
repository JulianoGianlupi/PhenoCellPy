import tissue_forge as tf
import numpy as np

import sys
from os.path import abspath

sys.path.extend([abspath("../")])  # todo: make this more refined

import Phenotypes as pheno


def get_radius_sphere(volume):
    return ((1 / (np.pi * 4 / 3)) * volume) ** (1 / 3)


# potential cutoff distance
cutoff = 1

# space set up
dim = [20, 20, 20]
tf.init(dim=dim)

pot = tf.Potential.morse(d=0.1, a=6, min=-1, max=1)

# Particle types

volume = 20

radius = get_radius_sphere(volume)

mass = volume

dt = 5  # min/time step

ki67_basic = pheno.cycles.Ki67Basic(dt=dt, target_volumes=[mass, mass], volumes=[mass, mass])


class CellType(tf.ParticleTypeSpec):
    mass = mass
    target_temperature = 0
    radius = radius
    pvolume = volume
    dynamics = tf.Overdamped
    cycle = ki67_basic


Cell = CellType.get()

tf.bind.types(pot, Cell, Cell)

first_cell = Cell([d // 2 for d in dim])
first_cell.cycle = ki67_basic


def step_cycle_and_divide(event):
    # print("Time step")
    for p in Cell.items():
        print("in particle loop")
        print(hasattr(p, "cycle"))
        p.cycle.current_phase.simulated_cell_volume = p.pvolume
        phase_change, death, division = p.cycle.time_step_cycle()
        print(phase_change, death, division)
        if phase_change and len(CellType.items()) < 10:
            print("@@@\nPHASE CHANGE\n@@@")

        radius = get_radius_sphere(p.cycle.current_phase.volume)

        # book-keeping, making sure the simulated cell grows
        if p.radius < radius:
            print("update volume")
            p.radius = radius
            p.mass = p.cycle.current_phase.volume
            p.volume = p.cycle.current_phase.volume

        # if division occurs, divide
        if division:
            # save cell attribs to halve later
            radius = p.radius
            mass = p.mass
            volume = p.pvolume

            # divide and reasign attribs (is this step necessary?)
            child = p.split()
            child.radius = p.radius = radius / 2
            child.mass = p.mass = mass / 2
            child.pvolume = p.pvolume = volume / 2
    return 0


tf.event.on_time(invoke_method=step_cycle_and_divide, period=tf.Universe.dt)

# run the simulator interactive
tf.run()
