
from numpy import exp

from numpy.random import uniform


class Phase:

    def __init__(self, index: int = None, previous_phase_index: int = None, next_phase_index: int = None,
                 time_unit: str = "min", name: str = None, division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 10,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None):

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
            self.transition_to_next_phase = self._transition_to_next_phase
        else:
            self.transition_to_next_phase = transition_to_next_phase

    def _transition_to_next_phase(self, dt):
        """
        Default stochastic phase transition function. Calculates a Poisson probability based on dt and
        self.phase_duration, rolls a random number, and returns random number < probability
        :param dt: float. Time-step length in the time units of Phase
        :return: bool. random number < probability of transition
        """
        prob = 1 - exp(-dt/self.phase_duration)
        return uniform() < prob

    def time_step_phase(self, dt):
        """

        :param dt: float. Time-step length in the time units of Phase
        :return: tuple. First element of tuple: bool denoting if the cell moves to the next phase. Second element:
        denotes if the cell leaves the cell cycle and enters quiescence.
        """
        if dt <= 0:
            raise ValueError(f"'dt' must be greater than 0. Got {dt}.")

        self.time_in_phase += dt

        if self.arrest_function is not None:
            if self.arrest_function(self.arrest_function_args):
                return False, True

        if self.fixed_duration and self.time_in_phase > self.phase_duration:
            if self.exit_function is not None:
                self.exit_function(self.exit_function_args)
            return True, False
        elif not self.fixed_duration:
            transition = self.transition_to_next_phase(dt)
            return transition, False
        return False, False


if __name__ == '__main__':
    pass

    # testPhase = Phase(phase_duration=-1)
