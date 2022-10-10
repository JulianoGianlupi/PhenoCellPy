
import tissue_forge as tf

import sys
from os.path import abspath
sys.path.extend([abspath("../")])  # todo: make this more refined

import Phenotypes as pheno

# potential cutoff distance
cutoff = 1

# space set up
dim = [20, 20, 20]
tf.init(dim=dim)


pot = tf.Potential.morse(d=0.1, a=6, min=-1, max=1)

# Particle types

mass = 20

dt = 5  # 5 min/time step

ki67_basic = pheno.cycles.Ki67Basic(dt=dt, target_volumes=[mass, mass], volumes=[mass, mass])

class CellType(tf.ParticleTypeSpec):
    mass = 20
    target_temperature = 0
    radius = 0.5
    dynamics = tf.Overdamped
    cycle = ki67_basic
    @staticmethod
    def on_register(ptype):
        def step_cell_cycle(event: tf.event.ParticleTimeEvent):
            m = event.targetParticle
            m.cycle.current_phase.simulated_cell_volume = m.mass
            phase_change, death, division = m.cycle.time_step_cycle()

            # print(f"Time-step cycle, {m.cycle.time_in_cycle}, {m.cycle.current_phase.time_in_phase}, "
            #       f"{m.cycle.current_phase.index}")
            print("Time step")

            if m.mass < m.cycle.current_phase.volume:
                m.mass = m.cycle.current_phase.volume

            if division:
                d = m.fission()

                m.radius = d.radius = CellType.radius
                m.mass = d.mass = CellType.mass

                m.cycle.current_phase.target_volume = d.cycle.current_phase.target_volume = CellType.mass
                m.cycle.current_phase.volume = d.cycle.current_phase.volume = CellType.mass
                m.cycle.current_phase.simulated_cell_volume = d.cycle.current_phase.simulated_cell_volume = \
                    CellType.mass
                print('fission:', len(event.targetType.items()))
        tf.event.on_particletime(ptype=ptype, invoke_method=step_cell_cycle, period=1)  # distribution=???


Cell = CellType.get()

tf.bind.types(pot, Cell, Cell)

Cell([d//2 for d in dim])

# run the simulator interactive
tf.run()

