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

necrosis_phenotype = pheno.phenotypes.get_phenotype_by_name("Standard necrosis model")

nuclear_initial_volume = .75 * mass
cytoplasm_initial_volume = .25 * mass

necrosis_phenotype = necrosis_phenotype(dt=dt,
                                        # the default volume
                                        # model uses a .75 fluid fraction
                                        nuclear_fluid=[.75 * nuclear_initial_volume, .75 * nuclear_initial_volume],
                                        nuclear_solid=[.25 * nuclear_initial_volume, .25 * nuclear_initial_volume],
                                        cytoplasm_fluid=[.75 * cytoplasm_initial_volume,
                                                         .75 * cytoplasm_initial_volume],
                                        cytoplasm_solid=[.25 * cytoplasm_initial_volume,
                                                         .25 * cytoplasm_initial_volume],
                                        simulated_cell_volume=mass)


class CellType(tf.ParticleTypeSpec):
    mass = mass
    target_temperature = 0
    radius = radius
    # volume = volume
    dynamics = tf.Overdamped


Cell = CellType.get()

tf.bind.types(pot, Cell, Cell)

rforce = tf.Force.random(mean=0, std=50)

# bind it just like any other force
tf.bind.force(rforce, Cell)
