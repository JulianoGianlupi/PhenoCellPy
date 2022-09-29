from numpy import exp

from numpy.random import uniform


# todo:
#  - finish generic cycle class
#  - implement phase transition
#  - implement quiescent phenotype
#  - implement physicell's phenotypes

class Phase:

    def __init__(self, index: int = None, previous_phase_index: int = None, next_phase_index: int = None,
                 dt: float = None, time_unit: str = "min", name: str = None, division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 10,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None):

        """
        :param index:
        :param previous_phase_index:
        :param next_phase_index:
        :param time_unit:
        :param name:
        :param division_at_phase_exit:
        :param removal_at_phase_exit:
        :param fixed_duration:
        :param phase_duration:
        :param entry_function:
        :param entry_function_args:
        :param exit_function:
        :param exit_function_args:
        :param arrest_function:
        :param arrest_function_args:
        :param transition_to_next_phase:
        """

        if index is None:
            self.index = 0  # int
        else:
            self.index = index

        self.previous_phase_index = previous_phase_index

        self.next_phase_index = next_phase_index

        self.time_unit = time_unit

        if dt <= 0 or dt is None:
            raise ValueError(f"'dt' must be greater than 0. Got {dt}.")
        self.dt = dt

        if name is None:
            self.name = "unnamed"  # string: phase's name
        else:
            self.name = name

        self.division_at_phase_exit = division_at_phase_exit  # bool flagging for division
        self.removal_at_phase_exit = removal_at_phase_exit  # bool flagging for removal (death?)

        self.fixed_duration = fixed_duration

        if phase_duration <= 0:
            raise ValueError(f"'phase_duration' must be greater than 0. Got {phase_duration}")

        self.phase_duration = phase_duration

        self.time_in_phase = 0

        self.entry_function = entry_function  # function to be executed upon entering this phase
        self.entry_function_args = entry_function_args

        self.exit_function = exit_function  # function to be executed just before exiting this phase
        self.exit_function_args = exit_function_args

        self.arrest_function = arrest_function  # function determining if cell will exit cell cycle and become quiescent
        self.arrest_function_args = arrest_function_args

        if transition_to_next_phase is None:
            self.custom_transition_function = False
            self.transition_to_next_phase_args = None
            if fixed_duration:
                self.transition_to_next_phase = self._transition_to_next_phase_deterministic
            else:
                self.transition_to_next_phase = self._transition_to_next_phase_stochastic
        else:
            self.custom_transition_function = True
            if type(transition_to_next_phase_args) != list:
                raise TypeError("Custom exit function selected but no args given. Was expecting "
                                f"'transition_to_next_phase_args' to be a list, got {transition_to_next_phase_args}.")
            self.transition_to_next_phase_args = transition_to_next_phase_args
            self.transition_to_next_phase = transition_to_next_phase

    def _transition_to_next_phase_stochastic(self):
        """
        Default stochastic phase transition function. Calculates a Poisson probability based on dt and
        self.phase_duration, rolls a random number, and returns random number < probability
        :param args: float. Time-step length in the time units of Phase
        :return: bool. random number < probability of transition
        """
        prob = 1 - exp(-self.dt / self.phase_duration)
        return uniform() < prob

    def _transition_to_next_phase_deterministic(self):
        return self.time_in_phase > self.phase_duration

    def time_step_phase(self):
        """

        :return: tuple. First element of tuple: bool denoting if the cell moves to the next phase. Second element:
        denotes if the cell leaves the cell cycle and enters quiescence.
        """
        self.time_in_phase += self.dt

        if self.arrest_function is not None:
            if self.arrest_function(self.arrest_function_args):
                return False, True

        if self.custom_transition_function:
            transition = self.transition_to_next_phase(self.transition_to_next_phase_args)
        else:
            transition = self.transition_to_next_phase()

        if transition and self.exit_function is not None:
            quies = self.exit_function(self.exit_function_args)
            return transition, quies
        elif transition:
            return transition, False

        return False, False


class QuiescentPhase(Phase):
    def __init__(self, index: int = 999, next_phase_index: int = 0, time_unit: str = "min", ):
        # todo: this whole class
        super().__init__(index=index, next_phase_index=next_phase_index)
        return


class Cycle:
    def __init__(self):
        self.phases = [Phase(previous_phase_index=0, next_phase_index=0)]
        self.quiecent_phase = QuiescentPhase()
        self.current_phase = self.phases[0]
        self.time_in_cycle = 0

    def time_step_cycle(self, dt):

        if dt <= 0:
            raise ValueError(f"'dt' must be greater than 0. Got {dt}.")

        self.time_in_cycle += dt

        next_phase, quies = self.current_phase.time_step_phase()

        if next_phase:
            changed_phases, cell_dies, cell_divides = self.go_to_next_phase()
            return changed_phases, cell_dies, cell_divides
        elif quies:
            self.go_to_quiescence()
            return True, False, False

    def go_to_next_phase(self):
        divides = self.current_phase.division_at_phase_exit
        dies = self.current_phase.removal_at_phase_exit
        self.set_phase(self.current_phase.next_phase_index)
        return True, dies, divides

    def set_phase(self, idx):
        self.current_phase = self.phases[idx]
        self.current_phase.time_in_phase = 0

    def go_to_quiescence(self):
        self.current_phase = self.quiecent_phase
        self.current_phase.time_in_phase = 0


class SimpleLiveCycle(Cycle):
    def __init__(self, time_unit: str = "min", name: str = "simple_live"):
        super().__init__()
        self.time_unit = time_unit

        self.name = name

        self.phases = [Phase(index=0, previous_phase_index=0, next_phase_index=0, time_unit=time_unit, name="alive",
                             division_at_phase_exit=True, )]


if __name__ == '__main__':
    pass

    # testPhase = Phase(phase_duration=-1)
