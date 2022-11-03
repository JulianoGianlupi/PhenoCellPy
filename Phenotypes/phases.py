from numpy import exp
from numpy.random import uniform

from Phenotypes.cell_volume import CellVolumes


# todo:
#  - keep updating the old volume to new volume

class Phase:
    """
    Base class to define phases of a cell cycle.

    Methods:
    -------

    time_step_phase()
        Time steps the phase. Returns a tuple (did the cell transition to the next phase, did the cell enter
        quiescence). See time_step_phase's documentation for further explanation.

    transition_to_next_phase() or transition_to_next_phase(*args)
        One of the default transition functions (`_transition_to_next_phase_deterministic`,
        `_transition_to_next_phase_stochastic`) or an user defined function. If user defined, it will be called with
        `transition_to_next_phase_args` as args. Must return a bool denoting if the transition occurs or not.

    _transition_to_next_phase_deterministic()
        Default deterministic transition function. Returns `time_in_phase > phase_duration`

    _transition_to_next_phase_stochastic()
        Default stochastic transition function. Probability of transition depends on `dt` and `phase_duration`

    entry_function(*args)
        Optional user defined function to be executed upon entering this phase. It gets called using attribute
        `entry_function_args`. Must have no return

    exit_function(*args)
        Optional user defined function to be executed upon exiting this phase. It gets called using attribute
        `exit_function_args`. Must have no return

    arrest_function(*args)
        Optional user defined function that returns true if the cell should exit the cell cycle and enter quiescence

    update_volume() or update_volume(*args)
        Function to update the volume of the cell. Can be user defined, in that case it has to use args and cannot have
        a return value

    Attributes
    ----------

    name : str
        Name of the phase

    index : int
        The index of the phase in the cycle (0 idexed)

    previous_phase_index : int
        Index of the phase preceding this one (usually index-1)

    next_phase_index : int
        Index of the phase preceding this one (usually index+1)

    time_unit : str
        Time unit used by the model. TODO: defines unit conversions

    dt : float
        Size of the time-step in units of `time_unit`

    division_at_phase_exit : bool
        Flag for cell division at the end of the phase

    removal_at_phase_exit : bool
        Flag for cell removal (death) at the end of the phase

    phase_duration : float
        Average amount of time the phase lasts (in units of `time_unit`)

    fixed_duration : bool
        Flag for determining if the phase ends after a set amount of time or if it progresses to the next phase
        in a stochastic manner with a rate `1/phase_duration`

    time_in_phase : float
        Time spent in this phase

    simulated_cell_volume : float
        The volume of the simulated cell object in this phase



    entry_function_args : list
        args for `entry_function`

    exit_function_args : list
        args for `exit_function`

    arrest_function_args : list
        args for `arrest_function`

    transition_to_next_phase_args : list
        args for `transition_to_next_phase` if a custom one is defined

    """

    def __init__(self, index: int = None, previous_phase_index: int = None, next_phase_index: int = None,
                 dt: float = None, time_unit: str = "min", name: str = None, division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 10,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=None,
                 nuclear_biomass_change_rate=None, calcification_rate=None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None):
        # todo: replace the old volume things with the new volume things
        """
        :param fluid_change_rate:
        :param target_fluid_fraction:
        :param nuclear_fluid:
        :param nuclear_solid:
        :param nuclear_solid_target:
        :param cytoplasm_fluid:
        :param cytoplasm_solid:
        :param cytoplasm_solid_target:
        :param target_cytoplasm_to_nuclear_ratio:
        :param calcified_fraction:
        :param calcified_fraction:
        :param calcification_rate:
        :param cytoplasm_biomass_change_rate:
        :param nuclear_biomass_change_rate:
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
            self.transition_to_next_phase_args = [None]
            if fixed_duration:
                self.transition_to_next_phase = self._transition_to_next_phase_deterministic
            else:
                self.transition_to_next_phase = self._transition_to_next_phase_stochastic
        else:
            if type(transition_to_next_phase_args) != list:
                raise TypeError("Custom exit function selected but no args given. Was expecting "
                                "'transition_to_next_phase_args' to be a list, got "
                                f"{type(transition_to_next_phase_args)}.")
            self.transition_to_next_phase_args = transition_to_next_phase_args
            self.transition_to_next_phase = transition_to_next_phase

        if simulated_cell_volume is None:
            self.simulated_cell_volume = 1
        else:
            self.simulated_cell_volume = simulated_cell_volume

        # the default rates are reference values for MCF-7, in 1/min
        if cytoplasm_biomass_change_rate is None:
            self.cytoplasm_biomass_change_rate = 0.27 / 60.0
        else:
            self.cytoplasm_biomass_change_rate = cytoplasm_biomass_change_rate

        if nuclear_biomass_change_rate is None:
            self.nuclear_biomass_change_rate = 0.33 / 60.0
        else:
            self.nuclear_biomass_change_rate = nuclear_biomass_change_rate
        if calcification_rate is None:
            self.calcification_rate = 0
        else:
            if calcification_rate < 0:
                raise ValueError(f"`calcification_rate` must be >= 0, got {calcification_rate}")
            self.calcification_rate = calcification_rate

        if fluid_change_rate is None:
            self.fluid_change_rate = 3.0 / 60.0
        else:
            self.fluid_change_rate = fluid_change_rate

        # self.volume = CellVolumes(cytoplasm=cytoplasm_volume, target_cytoplasm=cytoplasm_target_volume,
        #                           target_cytoplasm_fluid_fraction=cytoplasm_target_fluid_fraction,
        #                           nuclear=nuclear_volume, target_nuclear=nuclear_target_volume,
        #                           target_nuclear_fluid_fraction=nuclear_target_fluid_fraction,
        #                           calcified_fraction=calcified_fraction)

        self.volume = CellVolumes(target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                                  nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                                  cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                                  cytoplasm_solid_target=cytoplasm_solid_target,
                                  target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                                  calcified_fraction=calcified_fraction)

    def update_volume(self):
        self.volume.update_volume(self.dt, self.fluid_change_rate, self.nuclear_biomass_change_rate,
                                  self.cytoplasm_biomass_change_rate, self.calcification_rate)

    def _transition_to_next_phase_stochastic(self, none):
        """
        Default stochastic phase transition function.

        Calculates a Poisson probability based on dt and self.phase_duration (p=1-exp(-dt/phase_duration), rolls a
        random number, and returns random number < probability.

        :param none: Not used. Place holder in case of user defined function with args
        :return: bool. random number < probability of transition
        """

        # the approximation 1-exp(-x) ~ x can be used. That approximation has a difference of 0.005 at x=0.1, which I'd
        # find acceptable. TODO: implement a check on self.dt / self.phase_duration, if it is < .1 use the approximation

        prob = 1 - exp(-self.dt / self.phase_duration)
        return uniform() < prob

    def _transition_to_next_phase_deterministic(self, none):
        """
        Default deterministic phase transition function.

        If the time spent in this phase is greater than the phase duration, go to the next phase.

        :param none: Not used. Place holder in case of user defined function with args
        :return:
        """
        return self.time_in_phase > self.phase_duration

    def time_step_phase(self):
        """

        Time steps the phase.

        This function increments the `time_in_phase` by `dt`. Updates the cell volume. Checks if the cell arrests its
        cycle (i.e., leaves the cycle; goes to quiescence). If the cell doesn't quies, this function checks if the cell
        should transition to the next phase.

        :return: tuple. First element of tuple: bool denoting if the cell moves to the next phase. Second element:
        denotes if the cell leaves the cell cycle and enters quiescence.
        """
        self.time_in_phase += self.dt

        self.update_volume()

        if self.arrest_function is not None:
            if self.arrest_function(*self.arrest_function_args):
                return False, True

        transition = self.transition_to_next_phase(*self.transition_to_next_phase_args)

        if transition and self.exit_function is not None:
            quies = self.exit_function(*self.exit_function_args)
            return transition, quies
        elif transition:
            return transition, False

        return False, False

    def _double_target_volume(self, *none):
        self.volume.nuclear_solid_target *= 2
        self.volume.cytoplasm_solid_target *= 2

    def __str__(self):
        return f"{self.name} phase"


class QuiescentPhase(Phase):
    """Default Quiescent Phase. Inherits Phase()"""

    def __init__(self, index: int = 9999, previous_phase_index: int = None, next_phase_index: int = 0, dt: float = None,
                 time_unit: str = "min", name: str = "quiescent", division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 4.59 * 60,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=0, nuclear_biomass_change_rate=0,
                 calcification_rate=0, target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None,
                 nuclear_solid_target=None, cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate,
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate)
        return


class Ki67Negative(Phase):
    """

    Defines Ki 67- phase.

    TODO: More description here

    """

    def __init__(self, index: int = 0, previous_phase_index: int = 1, next_phase_index: int = 1, dt: float = 0.1,
                 time_unit: str = "min", name: str = "Ki 67-", division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 4.59 * 60,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=None,
                 nuclear_biomass_change_rate=None, calcification_rate=None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate,
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate)


class Ki67Positive(Phase):
    """

    Defines the simple Ki 67+ phase. Inherits Phase().

    Methods
    -------

    _standard_Ki67_entry_function
        Default standard entry function to this phase.

    """

    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 0, dt: float = 0.1,
                 time_unit: str = "min", name: str = "Ki 67+", division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = True, phase_duration: float = 15.5 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=None,
                 nuclear_biomass_change_rate=None, calcification_rate=None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None):

        if entry_function is None:
            entry_function = self._double_target_volume
            entry_function_args = [None]
        elif type(entry_function_args) != list:
            raise TypeError("'entry_function' was defined but no value for 'entry_function_args' was given. Expected "
                            f"list got {type(entry_function_args)}")

        if target_fluid_fraction > 1:
            target_fluid_fraction = 1
        elif target_fluid_fraction < 0:
            target_fluid_fraction = 0

        if cytoplasm_biomass_change_rate is None and cytoplasm_fluid is not None and cytoplasm_solid is not None:
            cytoplasm_biomass_change_rate = (cytoplasm_fluid + cytoplasm_solid) / (phase_duration / dt)

        elif cytoplasm_biomass_change_rate is None and cytoplasm_fluid is not None:
            cytoplasm_biomass_change_rate = cytoplasm_fluid / (phase_duration / dt)

        elif cytoplasm_biomass_change_rate is None and cytoplasm_solid is not None:
            cytoplasm_biomass_change_rate = cytoplasm_solid / (phase_duration / dt)

        elif cytoplasm_biomass_change_rate is None:
            cytoplasm_biomass_change_rate = 1
        else:
            cytoplasm_biomass_change_rate = cytoplasm_biomass_change_rate

        if nuclear_biomass_change_rate is None and cytoplasm_fluid is not None and cytoplasm_solid is not None:
            nuclear_biomass_change_rate = (nuclear_fluid + nuclear_solid) / (phase_duration / dt)

        elif nuclear_biomass_change_rate is None and cytoplasm_fluid is not None:
            nuclear_biomass_change_rate = nuclear_fluid / (phase_duration / dt)

        elif nuclear_biomass_change_rate is None and cytoplasm_solid is not None:
            nuclear_biomass_change_rate = nuclear_solid / (phase_duration / dt)

        elif nuclear_biomass_change_rate is None:
            nuclear_biomass_change_rate = 1
        else:
            nuclear_biomass_change_rate = nuclear_biomass_change_rate

        if fluid_change_rate is None and cytoplasm_fluid is not None and nuclear_fluid is not None:
            fluid_change_rate = (cytoplasm_fluid + nuclear_fluid) / (phase_duration / dt)
        elif fluid_change_rate is None and cytoplasm_fluid is not None:
            fluid_change_rate = cytoplasm_fluid / (phase_duration / dt)
        elif fluid_change_rate is None and nuclear_fluid is not None:
            fluid_change_rate = nuclear_fluid / (phase_duration / dt)
        else:
            fluid_change_rate = 1

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate,
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate)


class Ki67PositivePreMitotic(Ki67Positive):

    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 2, dt: float = 0.1,
                 time_unit: str = "min", name: str = "Ki 67+ pre-mitotic", division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = True, phase_duration: float = 13.0 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=None,
                 nuclear_biomass_change_rate=None, calcification_rate=None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate,
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate)


class Ki67PositivePostMitotic(Phase):
    def __init__(self, index: int = 2, previous_phase_index: int = 1, next_phase_index: int = 0, dt: float = 0.1,
                 time_unit: str = "min", name: str = "Ki 67+ post-mitotic", division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = True, phase_duration: float = 2.5 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=None,
                 nuclear_biomass_change_rate=None, calcification_rate=None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None):

        if entry_function is None:
            entry_function = self._standard_Ki67_positive_postmit_entry_function
            entry_function_args = [None]
        elif type(entry_function_args) != list:
            raise TypeError("'entry_function' was defined but no value for 'entry_function_args' was given. Expected "
                            f"list got {type(entry_function_args)}")

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate,
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate)

    def _standard_Ki67_positive_postmit_entry_function(self, *args):
        self.volume.cytoplasm_solid_target /= 2
        self.volume.nuclear_solid_target /= 2


class G0G1(Phase):
    def __init__(self, index: int = 0, previous_phase_index: int = 2, next_phase_index: int = 1, dt: float = 0.1,
                 time_unit: str = "min", name: str = "G0/G1", division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 5.15 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=None,
                 nuclear_biomass_change_rate=None, calcification_rate=None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate,
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate)


class S(Phase):
    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 2, dt: float = 0.1,
                 time_unit: str = "min", name: str = "S", division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 8 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=None,
                 nuclear_biomass_change_rate=None, calcification_rate=None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None):

        if entry_function is None:
            entry_function = self._double_target_volume
            entry_function_args = [None]
        elif type(entry_function_args) != list:
            raise TypeError("'entry_function' was defined but no value for 'entry_function_args' was given. Expected "
                            f"list got {type(entry_function_args)}")

        if cytoplasm_biomass_change_rate is None and cytoplasm_fluid is not None and cytoplasm_solid is not None:
            cytoplasm_biomass_change_rate = (cytoplasm_fluid + cytoplasm_solid) / (phase_duration / dt)

        elif cytoplasm_biomass_change_rate is None and cytoplasm_fluid is not None:
            cytoplasm_biomass_change_rate = cytoplasm_fluid / (phase_duration / dt)

        elif cytoplasm_biomass_change_rate is None and cytoplasm_solid is not None:
            cytoplasm_biomass_change_rate = cytoplasm_solid / (phase_duration / dt)

        elif cytoplasm_biomass_change_rate is None:
            cytoplasm_biomass_change_rate = 1
        else:
            cytoplasm_biomass_change_rate = cytoplasm_biomass_change_rate

        if nuclear_biomass_change_rate is None and cytoplasm_fluid is not None and cytoplasm_solid is not None:
            nuclear_biomass_change_rate = (nuclear_fluid + nuclear_solid) / (phase_duration / dt)

        elif nuclear_biomass_change_rate is None and cytoplasm_fluid is not None:
            nuclear_biomass_change_rate = nuclear_fluid / (phase_duration / dt)

        elif nuclear_biomass_change_rate is None and cytoplasm_solid is not None:
            nuclear_biomass_change_rate = nuclear_solid / (phase_duration / dt)

        elif nuclear_biomass_change_rate is None:
            nuclear_biomass_change_rate = 1
        else:
            nuclear_biomass_change_rate = nuclear_biomass_change_rate

        if fluid_change_rate is None and cytoplasm_fluid is not None and nuclear_fluid is not None:
            fluid_change_rate = (cytoplasm_fluid + nuclear_fluid) / (phase_duration / dt)
        elif fluid_change_rate is None and cytoplasm_fluid is not None:
            fluid_change_rate = cytoplasm_fluid / (phase_duration / dt)
        elif fluid_change_rate is None and nuclear_fluid is not None:
            fluid_change_rate = nuclear_fluid / (phase_duration / dt)
        else:
            fluid_change_rate = 1

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate,
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate)


class G2M(Phase):
    def __init__(self, index: int = 2, previous_phase_index: int = 1, next_phase_index: int = 0, dt: float = 0.1,
                 time_unit: str = "min", name: str = "G2/M", division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 5 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=None,
                 nuclear_biomass_change_rate=None, calcification_rate=None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate,
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate)


class Apoptosis(Phase):
    def __init__(self, index: int = 0, previous_phase_index: int = 0, next_phase_index: int = 0, dt: float = 0.1,
                 time_unit: str = "min", name: str = "Apoptosis", division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = True, fixed_duration: bool = True, phase_duration: float = 8.6 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate: float = 1 / 60,
                 nuclear_biomass_change_rate: float = 0.35 / 60, unlysed_fluid_change_rate: float = 3 / 60,
                 lysed_fluid_change_rate: float = 0, calcification_rate: float = 0, relative_rupture_volume: float = 2,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None):
        # todo: figure out what unlysed_fluid_change_rate, lysed_fluid_change_rate, relative_rupture_volume are supposed
        #  to do

        if entry_function is None:
            entry_function = self._standard_apoptosis_entry

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate,
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate)

    def _standard_apoptosis_entry(self):
        self.volume.target_cytoplasm = 0
        self.volume.target_nuclear = 0


if __name__ == '__main__':
    test_ki = Ki67Positive(dt=0.1)
    print(test_ki.index)
