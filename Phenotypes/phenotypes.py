# todo:
#  - finish generic cycle class
#  - implement phase transition
#  - implement quiescent phenotype
#  - implement physicell's phenotypes
#  - documentation
#  - functions to attach to object

import Phenotypes.cycles as cycles

if __name__ == '__main__':
    pass

    testCycle = cycles.Ki67Basic()
    print(testCycle.time_in_cycle, testCycle.current_phase.name, testCycle.current_phase.time_in_phase)
    for i in range(10000):
        changed_phase, died, divides = testCycle.time_step_cycle()

        if changed_phase or died or divides:
            print(changed_phase, died, divides)
            print(testCycle.time_in_cycle, testCycle.current_phase.name, testCycle.current_phase.time_in_phase)

    pass
