import tissue_forge as tf
import numpy as np

import sys
from os.path import abspath
# import keyboard
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

# pot = tf.Potential.morse(d=0.1, a=6, min=-1, max=1)
pot = tf.Potential.glj(1)

# Particle types

# volume = 10

mass = 40

radius = .4 #get_radius_sphere(volume)

global density
density = mass / ((4/3)*np.pi*radius*radius*radius)



dt = 1  # min/time step

ki67_basic = pheno.cycles.Ki67Basic(dt=dt, target_volumes=[mass, mass], volumes=[mass, mass])


class CellType(tf.ParticleTypeSpec):
    mass = mass
    target_temperature = 0
    radius = radius
    # volume = volume
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
    # print("Time step")
    for p in Cell.items():
        print(p.id)
        # print("in particle loop")
        # print(hasattr(p, "volume"))
        pcycle = cells_cycles[f"{p.id}"]
        pcycle.current_phase.simulated_cell_volume = p.mass * density
        phase_change, death, division = pcycle.time_step_cycle()
        # print(pcycle.current_phase)
        # print(pcycle.current_phase, phase_change, death, division, pcycle.current_phase.time_in_phase)
        # print(1 - np.exp(-pcycle.current_phase.dt / pcycle.current_phase.phase_duration))
        if phase_change and len(Cell.items()) < 10:
            print("@@@\nPHASE CHANGE\n@@@")
            # time.sleep(1)

        radius = get_radius_sphere(pcycle.current_phase.volume)

        # book-keeping, making sure the simulated cell grows
        if p.radius < radius:
            # print("update volume")
            p.radius = radius
            p.mass = ((4/3)*np.pi*radius*radius*radius) * density
            # p.volume = p.cycle.current_phase.volume

        # if division occurs, divide
        if division:
            # save cell attribs to halve later
            # radius = p.radius
            cur_mass = p.mass
            # volume = p.pvolume

            # divide and reasign attribs (is this step necessary?)
            child = p.split()

            cells_cycles[f"{child.id}"] = ki67_basic

            child.mass = p.mass = cur_mass / 2
            child.radius = p.radius = get_radius_sphere((cur_mass / 2)/density)

            cells_cycles[f"{child.id}"].volume = child.mass*density
            cells_cycles[f"{child.id}"].simulated_cell_volume = child.mass*density

            # child.pvolume = p.pvolume = volume / 2
    return 0


tf.event.on_time(invoke_method=step_cycle_and_divide, period=tf.Universe.dt)

# run the simulator interactive
tf.run()
