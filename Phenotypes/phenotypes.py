# todo:
#  - finish generic cycle class
#  - implement phase transition
#  - implement quiescent phenotype
#  - implement physicell's phenotypes
#  - documentation
#  - functions to attach to object

from Phenotypes.cycles import SimpleLiveCycle

if __name__ == '__main__':
    pass

    testCycle = SimpleLiveCycle()
    for i in range(1000):
        changed_phase, died, divides = testCycle.time_step_cycle()
        print(testCycle.time_in_cycle, testCycle.current_phase.name, testCycle.current_phase.time_in_phase)
        if changed_phase or died or divides:
            print(changed_phase, died, divides)
    pass
