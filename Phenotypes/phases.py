from numpy import exp
from numpy.random import uniform


class Phase:

    def __init__(self, index: int = None, previous_phase_index: int = None, next_phase_index: int = None,
                 dt: float = None, time_unit: str = "min", name: str = None, division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 10,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None, target_volume: float = None,
                 volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None):

        """
        :param update_volume_rate:
        :param target_volume:
        :param volume:
        :param update_volume:
        :param update_volume_args:
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

        if dt is None or dt <= 0:
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
        if self.exit_function is not None and type(self.exit_function_args) != list:
            raise TypeError("Exit function defined but no args given. Was expecting "
                            f"'exit_function_args' to be a list, got {type(exit_function_args)}.")

        self.arrest_function = arrest_function  # function determining if cell will exit cell cycle and become quiescent
        self.arrest_function_args = arrest_function_args

        if self.arrest_function is not None and type(self.arrest_function_args) != list:
            raise TypeError("Arrest function defined but no args given. Was expecting "
                            f"'arrest_function_args' to be a list, got {type(arrest_function_args)}.")

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
                                "'transition_to_next_phase_args' to be a list, got "
                                f"{type(transition_to_next_phase_args)}.")
            self.transition_to_next_phase_args = transition_to_next_phase_args
            self.transition_to_next_phase = transition_to_next_phase

        if volume is None:
            self.volume = 1
        else:
            self.volume = volume

        if target_volume is None:
            self.target_volume = 1
        else:
            self.target_volume = target_volume

        if update_volume_rate is None:
            self.update_volume_rate = 0.1
        else:
            self.update_volume_rate = update_volume_rate

        if update_volume is None:
            self.update_volume = self._update_volume
            self.update_volume_args = None
        else:
            self.update_volume = update_volume
            self.update_volume_args = update_volume_args

    def _update_volume(self, none):
        self.volume += self.update_volume_rate * (self.target_volume - self.volume)

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

        if self.update_volume:
            self.update_volume(self.update_volume_args)

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
    def __init__(self, index: int = 9999, next_phase_index: int = 0, time_unit: str = "min", dt: float = None,
                 fixed_duration: bool = True, phase_duration: float = 4.59 * 60, transition_to_next_phase=None,
                 transition_to_next_phase_args: list = None, exit_function=None,
                 exit_function_args: list = None, update_volume=False, volume=None, target_volume=None):
        super().__init__(index=index, next_phase_index=next_phase_index, time_unit=time_unit, dt=dt,
                         fixed_duration=fixed_duration, phase_duration=phase_duration,
                         transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, update_volume=update_volume, volume=volume,
                         target_volume=target_volume)
        return


class Ki67Negative(Phase):
    def __init__(self, name="Ki 67 negative", dt=0.1, time_unit="min", phase_duration=4.59 * 60, fixed_duration=False,
                 index=0, next_phase_index=1, previous_phase_index=1):
        super().__init__(name=name, dt=dt, time_unit=time_unit, phase_duration=phase_duration,
                         fixed_duration=fixed_duration, index=index, next_phase_index=next_phase_index,
                         previous_phase_index=previous_phase_index)


class Ki67Positive(Phase):
    def __init__(self, index=None, previous_phase_index=None, next_phase_index=None, dt=None, time_unit="min",
                 name=None, division_at_phase_exit=True, removal_at_phase_exit=False, fixed_duration=False,
                 entry_function=None, entry_function_args=None, phase_duration=10):
        if entry_function is None:
            entry_function = self._standard_Ki67_entry_function
            entry_function_args = [None]
        elif type(entry_function_args) != list:
            raise TypeError("'entry_function' was defined but no value for 'entry_function_args' was given. Expected "
                            f"list got {type(entry_function_args)}")

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit)

    def _standard_Ki67_entry_function(self, *args):
        self.target_volume *= 2


if __name__ == '__main__':
    test_ki = Ki67Positive(dt=0.1)
    print(test_ki.index)