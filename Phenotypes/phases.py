from numpy import exp
from numpy.random import uniform
from warnings import warn


class CellVolumes:
    # total
    total = 0

    # fluid volumes
    fluid = 0
    fluid_change_rate = 0
    target_fluid_fraction = 0
    fluid_fraction = 0

    # nuclear volumes
    nuclear = 0
    nuclear_solid = 0
    target_solid_nuclear = 0
    nuclear_biomass_change_rate = 0

    # cytoplasm volumes
    cytoplasm = 0
    cyto_nucl_ratio = 1
    target_solid_cytoplasm = 0
    cytoplasm_solid = 0
    cytoplasm_fluid = 0
    cytoplasm_biomass_change_rate = 0

    # calcified
    calcified_fraction = 0
    calcification_rate = 0


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

    volume : float
        The volume of the cell in this phase

    target_volume : float
        The target volume of the cell in this phase

    cytoplasmic_volume : float
        The cytoplasm's volume

    cytoplasmic_solid_fraction: float
        The ratio of cytoplasmic volume that should be solid, 0<=cytoplasmic_solid_fraction<=1

    cyto_nucl_ratio: float
        The volume of the cytoplasm in relation to the nuclear volume.
        cytoplasm volume =  cyto_nucl_ratio * nuclear_volume

    nuclear_solid_fraction: float
        The ratio of nucler volume that should be solid, 0<=nuclear_solid_fraction<=1

    simulated_cell_volume : float
        The volume of the simulated cell object in this phase

    update_volume_rate : float
        Amount of volume the cell will grow by if it is below the target and the default update volume function is used.

    update_volume_args : list
        Args used by custom `update_volume` function

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
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None, target_volume: float = None,
                 volume: float = None, cytoplasmic_ratio=None, cytoplasmic_solid_fraction=None, cyto_nucl_ratio=None,
                 nuclear_solid_fraction=None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None):

        """
        :param cytoplasmic_ratio:
        :param cytoplasmic_solid_fraction:
        :param cyto_nucl_ratio:
        :param nuclear_solid_fraction:
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

        self.new_volume = CellVolumes()

        if volume is None:
            self.volume = 1
        else:
            self.volume = volume

        if target_volume is None:
            self.target_volume = 1
        else:
            self.target_volume = target_volume

        if cytoplasmic_ratio is None:
            self.cytoplasm_ratio = 0.9
        else:
            self.cytoplasm_ratio = cytoplasmic_ratio

        if cytoplasmic_solid_fraction is None:
            self.cytoplasmic_solid_fraction = 0
        else:
            self.cytoplasmic_solid_fraction = cytoplasmic_solid_fraction

        if cyto_nucl_ratio is None:
            self.cyto_nucl_ratio = (1 / self.cytoplasm_ratio) - 1
        else:
            self.cyto_nucl_ratio = cyto_nucl_ratio

        if nuclear_solid_fraction is None:
            self.nuclear_solid_fraction = 0
        else:
            self.nuclear_solid_fraction = nuclear_solid_fraction

        if simulated_cell_volume is None:
            self.simulated_cell_volume = 1
        else:
            self.simulated_cell_volume = simulated_cell_volume

        if update_volume_rate is None:
            self.update_volume_rate = 1
        else:
            self.update_volume_rate = update_volume_rate

        if update_volume is None:
            self.update_volume = self._update_volume
            self.update_volume_args = [None]
        else:
            self.update_volume = update_volume
            self.update_volume_args = update_volume_args

    @property
    def cytoplasm_ratio(self):
        return self.__cyto_ratio

    @cytoplasm_ratio.setter
    def cytoplasm_ratio(self, value):
        if not 0 <= value <= 1:
            raise ValueError(f"`cytoplasmic_ratio` must be in the range [0, 1]. "
                             f"Got {value}")
        self.__cyto_ratio = value

    @property
    def cyto_nucl_ratio(self):
        return self.__cyto_nucl_ratio

    @cyto_nucl_ratio.setter
    def cyto_nucl_ratio(self, value):
        if value < 0:
            raise ValueError(f"`cyto_nucl_ratio` must be >=0")
        self.__cyto_nucl_ratio = value

    @property
    def cytoplasmic_solid_fraction(self):
        return self.__csf

    @cytoplasmic_solid_fraction.setter
    def cytoplasmic_solid_fraction(self, value):
        if not 0 <= value <= 1:
            raise ValueError(f"`cytoplasmic_solid_fraction` must be in the range [0, 1]. "
                             f"Got {value}")

        if self.cytoplasm_ratio * self.target_volume + \
                self.cyto_nucl_ratio * self.cytoplasm_ratio * self.target_volume > self.target_volume:
            message = "WARNING: `cytoplasmic_volume` + `cyto_nucl_ratio` * `cytoplasmic_volume` > `target_volume`. " \
                      "Resetting `target_volume` = `cytoplasmic_volume` + `cyto_nucl_ratio` * `cytoplasmic_volume`"
            warn(message)
            # target volume = target cyto volume + target nucl volume =
            #     = cytoplasm_ratio * target volume + cyto_nucl_ratio * cytoplasm_ratio * target volume =
            #     = target volume * (cytoplasm_ratio + cyto_nucl_ratio * cytoplasm_ratio)
            #     = target volume * (cytoplasm_ratio * (1 + cyto_nucl_ratio))
            self.target_volume = ((1 + self.cyto_nucl_ratio) * self.cytoplasm_ratio) * self.target_volume

        self.__csf = value

    @property
    def nuclear_solid_fraction(self):
        return self.__nsf

    @nuclear_solid_fraction.setter
    def nuclear_solid_fraction(self, value):
        if not 0 <= value <= 1:
            raise ValueError(f"`nuclear_solid_fraction` must be in the range [0,1]. Got {value}")
        self.__nsf = value

    @property
    def cytoplasm_volume(self):
        return self.cytoplasm_ratio * self.volume

    @property
    def cytoplasm_solid_volume(self):
        return self.cytoplasmic_solid_fraction * self.cytoplasm_volume

    @property
    def cytoplasm_fluid_volume(self):
        return (1 - self.cytoplasmic_solid_fraction) * self.cytoplasm_volume

    @property
    def nuclear_volume(self):
        return self.cyto_nucl_ratio * self.cytoplasm_volume

    @property
    def nuclear_solid_volume(self):
        return self.nuclear_solid_fraction * self.nuclear_volume

    @property
    def nuclear_fluid_volume(self):
        return (1 - self.nuclear_solid_fraction) * self.nuclear_volume

    def _update_volume(self, none):
        """
        Updates the volume if cell is below target.

        If the Phase's volume is below the Phase's target this function increments `self.volume` by
        `self.update_volume_rate`.

        :param none: Not used. Place holder in case of user defined function with args
        :return:
        """

        self.volume += self.dt * self.update_volume_rate * (self.target_volume - self.volume)

        if self.volume < 0:
            self.volume = 0

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

        if self.update_volume:
            self.update_volume(*self.update_volume_args)

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

    def _double_target_volume(self):
        self.target_volume *= 2

    def __str__(self):
        return f"{self.name} phase"


class QuiescentPhase(Phase):
    """Default Quiescent Phase. Inherits Phase()"""

    def __init__(self, index: int = 9999, previous_phase_index: int = None, next_phase_index: int = 0,
                 dt: float = None, time_unit: str = "min", name: str = "quiescent",
                 division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 4.59 * 60,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 target_volume: float = None, volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args, target_volume=target_volume,
                         volume=volume, update_volume=update_volume, update_volume_args=update_volume_args,
                         update_volume_rate=update_volume_rate, simulated_cell_volume=simulated_cell_volume)
        return


class Ki67Negative(Phase):
    """

    Defines Ki 67- phase.

    TODO: More description here

    """

    def __init__(self, index: int = 0, previous_phase_index: int = 1, next_phase_index: int = 1,
                 dt: float = 0.1, time_unit: str = "min", name: str = "Ki 67-",
                 division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 4.59 * 60,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 target_volume: float = None,
                 volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args, target_volume=target_volume,
                         volume=volume, update_volume=update_volume, update_volume_args=update_volume_args,
                         update_volume_rate=update_volume_rate, simulated_cell_volume=simulated_cell_volume)


class Ki67Positive(Phase):
    """

    Defines the simple Ki 67+ phase. Inherits Phase().

    Methods
    -------

    _standard_Ki67_entry_function
        Default standard entry function to this phase.

    """

    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 0,
                 dt: float = 0.1, time_unit: str = "min", name: str = "Ki 67+",
                 division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = True, phase_duration: float = 15.5 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None, target_volume: float = None,
                 volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None):

        if entry_function is None:
            entry_function = self._double_target_volume
            entry_function_args = [None]
        elif type(entry_function_args) != list:
            raise TypeError("'entry_function' was defined but no value for 'entry_function_args' was given. Expected "
                            f"list got {type(entry_function_args)}")

        if update_volume_rate is None:
            update_volume_rate = target_volume / (phase_duration / dt)

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args, target_volume=target_volume,
                         volume=volume, update_volume=update_volume, update_volume_args=update_volume_args,
                         update_volume_rate=update_volume_rate, simulated_cell_volume=simulated_cell_volume)


class Ki67PositivePreMitotic(Ki67Positive):

    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 2,
                 dt: float = 0.1, time_unit: str = "min", name: str = "Ki 67+ pre-mitotic",
                 division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = True, phase_duration: float = 13.0 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None, target_volume: float = None,
                 volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, target_volume=target_volume, volume=volume,
                         update_volume=update_volume, update_volume_args=update_volume_args,
                         update_volume_rate=update_volume_rate, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args,
                         simulated_cell_volume=simulated_cell_volume, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args)


class Ki67PositivePostMitotic(Phase):
    def __init__(self, index: int = 2, previous_phase_index: int = 1, next_phase_index: int = 0,
                 dt: float = 0.1, time_unit: str = "min", name: str = "Ki 67+ post-mitotic",
                 division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = True, phase_duration: float = 2.5 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None, target_volume: float = None,
                 volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None):

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
                         transition_to_next_phase_args=transition_to_next_phase_args, target_volume=target_volume,
                         volume=volume, update_volume=update_volume, update_volume_args=update_volume_args,
                         update_volume_rate=update_volume_rate, simulated_cell_volume=simulated_cell_volume)

    def _standard_Ki67_positive_postmit_entry_function(self, *args):
        self.target_volume /= 2


class G0G1(Phase):
    def __init__(self, index: int = 0, previous_phase_index: int = 2, next_phase_index: int = 1,
                 dt: float = 0.1, time_unit: str = "min", name: str = "G0/G1",
                 division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 5.15 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None, target_volume: float = None,
                 volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args, target_volume=target_volume,
                         volume=volume, update_volume=update_volume, update_volume_args=update_volume_args,
                         update_volume_rate=update_volume_rate, simulated_cell_volume=simulated_cell_volume)


class S(Phase):
    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 2,
                 dt: float = 0.1, time_unit: str = "min", name: str = "S",
                 division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 8 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None, target_volume: float = None,
                 volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None):

        if entry_function is None:
            entry_function = self._double_target_volume
            entry_function_args = [None]
        elif type(entry_function_args) != list:
            raise TypeError("'entry_function' was defined but no value for 'entry_function_args' was given. Expected "
                            f"list got {type(entry_function_args)}")
        if update_volume_rate is None:
            update_volume_rate = target_volume / (phase_duration / dt)

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args, target_volume=target_volume,
                         volume=volume, update_volume=update_volume, update_volume_args=update_volume_args,
                         update_volume_rate=update_volume_rate, simulated_cell_volume=simulated_cell_volume)


class G2M(Phase):
    def __init__(self, index: int = 2, previous_phase_index: int = 1, next_phase_index: int = 0,
                 dt: float = 0.1, time_unit: str = "min", name: str = "G2/M",
                 division_at_phase_exit: bool = True,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 5 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None, target_volume: float = None,
                 volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args, target_volume=target_volume,
                         volume=volume, update_volume=update_volume, update_volume_args=update_volume_args,
                         update_volume_rate=update_volume_rate, simulated_cell_volume=simulated_cell_volume)


class Apoptosis(Phase):
    def __init__(self, index: int = 0, previous_phase_index: int = 0, next_phase_index: int = 0,
                 dt: float = 0.1, time_unit: str = "min", name: str = "Apoptosis",
                 division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = True, fixed_duration: bool = True, phase_duration: float = 8.6 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None, target_volume: float = 0,
                 volume: float = None, update_volume=None, update_volume_args: list = None,
                 update_volume_rate: float = None, simulated_cell_volume: float = None,
                 cytoplasmic_biomass_change_rate: float = 1 / 60, nuclear_biomass_change_rate: float = 0.35 / 60,
                 unlysed_fluid_change_rate: float = 3 / 60, lysed_fluid_change_rate: float = 0,
                 calcification_rate: float = 0, relative_rupture_volume: float = 2):
        # todo: define change volume rate

        if entry_function is None:
            entry_function = self._standard_apoptosis_entry

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args, transition_to_next_phase=transition_to_next_phase,
                         transition_to_next_phase_args=transition_to_next_phase_args, target_volume=target_volume,
                         volume=volume, update_volume=update_volume, update_volume_args=update_volume_args,
                         update_volume_rate=update_volume_rate, simulated_cell_volume=simulated_cell_volume)

    def _standard_apoptosis_entry(self):
        return


if __name__ == '__main__':
    test_ki = Ki67Positive(dt=0.1)
    print(test_ki.index)
