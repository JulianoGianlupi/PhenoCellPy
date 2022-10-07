import Phenotypes.phases as Phases


class Cycle:
    """

    Base class to define a cell cycle.

    TODO: more description


    Methods:
    --------

    time_step_cycle()
        Time-steps the cycle model. Returns a tuple (cell changed phases, cell died, cell divided). See
        :func:`time_step_cycle` for further information.

    go_to_next_phase()
        Goes to the next phase in the cycle

    set_phase(index)
        Sets the current cycle's phase to be phase of index `index`

    go_to_quiescence()
        Moves cycle to quiescent phase

    Attributes
    ----------

    name : str
        Name of the cycle

    dt : float
        Time-step size (in units of `time_unit`). Must be >0.

    time_unit : str
        Time unit. TODO: Defines time conversions

    phases : list
        Ordered list of phases this cycle goes through. Must be a list of :class:`Phases.Phase` objects.

    quiescent_phase : :class:`Phases.Phase`, None, or False
        If false the cycle won't have a (stand-alone) quiescent phase defined. If None the default
        :class:`Phases.QuiescentPhase` will be used as the stand-alone quiescent phase. If it is a :class:`Phases.Phase`
        object it will be used as the stand-alone quiescent phase.

    current_phase : :class:`Phases.Phase`
        The current (active) phase of the cycle.

    time_in_cycle : float
        Total time elapsed for the cycle

    """

    def __init__(self, name: str = "unnamed", dt: float = 1, time_unit: str = "min", phases: list = None,
                 quiescent_phase: Phases.Phase or False = None):

        self.name = name

        self.time_unit = time_unit

        if dt <= 0 or dt is None:
            raise ValueError(f"'dt' must be greater than 0. Got {dt}.")
        self.dt = dt
        if phases is None:
            self.phases = [Phases.Phase(previous_phase_index=0, next_phase_index=0, dt=self.dt, time_unit=time_unit)]
        else:
            self.phases = phases
        if quiescent_phase is None:
            self.quiescent_phase = Phases.QuiescentPhase(dt=self.dt)
        elif quiescent_phase is not None and not quiescent_phase:
            self.quiescent_phase = False
        elif not isinstance(quiescent_phase, Phases.Phase):
            raise ValueError(f"`quiescent_phase` must Phases.Phase object, False, or None. Got {quiescent_phase}")
        else:
            self.quiescent_phase = quiescent_phase
        self.current_phase = self.phases[0]
        self.time_in_cycle = 0

    def time_step_cycle(self):
        """
        Time-steps the cycle.

        Increments :attr:`time_in_cycle` by :attr:`dt`. Calls :func:`current_phase.time_step_phase`. If the phase time-
        step determines the cycle moves to the next phase (i.e., returns `True` for `next_phase`), calls
        :func:go_to_next_phase. If the phase time-step determines the cell exits the cell cycle and goes to quiescence
        (i.e., returns `True` for `quies`) calls :func:`go_to_quiescence`.

        :return: Flags (bool) for phase changing, cell death, and cell division
        :rtype: tuple of bool
        """

        self.time_in_cycle += self.dt

        next_phase, quies = self.current_phase.time_step_phase()

        if next_phase:
            changed_phases, cell_dies, cell_divides = self.go_to_next_phase()
            return changed_phases, cell_dies, cell_divides
        elif quies:
            self.go_to_quiescence()
            return True, False, False
        return False, False, False

    def go_to_next_phase(self):
        """

        Sets the cycle phase to be the phase of index :attr:`current_phase.next_phase_index`.

        Gets the flags for division and death upon phase exit. Sets the cycle phase to be
        :attr:`current_phase.next_phase_index` by calling :func:`set_phase`. Calls :func:`current_phase.entry_function`
        (after phase change). Returns that phase change has occurred, and the division and death upon phase exit flags.

        :return: Flags (bool) for phase changing, cell death, and cell division
        :rtype: tuple of bool
        """
        divides = self.current_phase.division_at_phase_exit
        dies = self.current_phase.removal_at_phase_exit
        self.set_phase(self.current_phase.next_phase_index)
        if self.current_phase.entry_function is not None:
            self.current_phase.entry_function(self.current_phase.entry_function_args)
        # self.current_phase.time_in_phase = 0
        return True, dies, divides

    def set_phase(self, idx):
        """

        Sets cycle phase to be phase of index :param:`idx`.

        This function moves the cycle to the next phase and coordinates their volume attributes. Saves current
        :attr:`current_phase.volume` and :attr:`current_phase.volume` to variables, sets :attr:`current_phase` to be
        `phases[idx]`, resets :attr:`current_phase.volume` and :attr:`current_phase.volume` to be the previously saved
        values, sets :attr:`current_phase.time_in_phase` to 0.

        :param idx: index of list :attr:`phases`, which phase to go to.
        :return: No return
        """
        volume = self.current_phase.volume
        tg_volume = self.current_phase.target_volume
        self.current_phase = self.phases[idx]
        self.current_phase.volume = volume
        self.current_phase.target_volume = tg_volume
        self.current_phase.time_in_phase = 0

    def go_to_quiescence(self):
        """

        Sets cycle phase to be the quiescent phase.

        This function checks that :attr:`quiescent_phase` is a :class:`Phases.Phase`, if that's the case it sets
        :attr:`current_phase` to be :attr:`quiescent_phase`. It also sets (after phase change) the
        :attr:`current_phase.volume` and :attr:`current_phase.volume` to be the previous phase :attr:`volume` by first
        saving it to a temporary variable. It resets :attr:`current_phase.time_in_phase` to be 0.

        :return: No return
        """
        if not isinstance(self.quiescent_phase, Phases.Phase):
            return
        volume = self.current_phase.volume
        self.quiescent_phase.volume = volume
        self.quiescent_phase.target_volume = volume
        self.current_phase = self.quiescent_phase
        self.current_phase.time_in_phase = 0

    def __str__(self):
        phases = ""
        for p in self.phases:
            phases += f"{p}, "
        if len(phases) > 2:
            phases = phases[:-2]
        return f"{self.name} cycle, phases: {phases}"


class SimpleLiveCycle(Cycle):
    """
        Simplest life cycle, it has only one phase. When "progressing" to the next phase it divides.
        """

    def __init__(self, time_unit: str = "min", name: str = "Simple Live", dt=1):
        phases = [Phases.Phase(index=0, previous_phase_index=0, next_phase_index=0, dt=dt, time_unit=time_unit,
                               name="alive", division_at_phase_exit=True, phase_duration=60)]
        super().__init__(name=name, time_unit=time_unit, phases=phases, quiescent_phase=False, dt=dt)


class Ki67Basic(Cycle):
    """

    Simple proliferating-quiescent phase. Does not use the stand-alone quiescent phase. Cell divides upon leaving Ki67+

    TODO: find where physicell found the definition of this cycle. Find their explanation for this cycle
    This is a two phase cycle. Ki67- (defined in :class:`Phases.Ki67Negative`) is the quiescent phase, Ki67-'s mean du-
    ration is 4.59h (stochastic transition to Ki67+ by default). Ki67+ (defined in :class:`Phases.Ki67Positive`) is the
    proliferating phase. It is responsible for doubling the volume of the cell at a rate of
    [increase in volume]/[Ki67+ duration]. Ki67+ duration is fixed (by default) at 15.5 hours. Once the cell exits Ki67+
    it divides.

    """

    def __init__(self, name="Ki67 Basic", dt=0.1, time_unit="min", quiescent_phase=False,
                 division_at_phase_exits=(False, True), removal_at_phase_exits=(False, False),
                 fixed_durations=(False, True), phase_durations: list = (4.59 * 60, 15.5 * 60.0),
                 entry_functions=(None, None), entry_functions_args=(None, None), exit_functions=(None, None),
                 exit_functions_args=(None, None), arrest_functions=(None, None), arrest_functions_args=(None, None),
                 transitions_to_next_phase=(None, None), transitions_to_next_phase_args: list = (None, None),
                 target_volumes: list = (1, 1), volumes: list = (1, 1), update_volumes=(None, None),
                 update_volumes_args: list = (None, None), update_volume_rates=(None, None),
                 simulated_cell_volume=None):

        self._check_arguments(2, name, target_volumes, division_at_phase_exits, removal_at_phase_exits, fixed_durations,
                              phase_durations, entry_functions, entry_functions_args, exit_functions,
                              exit_functions_args, arrest_functions, arrest_functions_args, transitions_to_next_phase,
                              transitions_to_next_phase_args, update_volumes, update_volumes_args, update_volume_rates)

        Ki67_positive = Phases.Ki67Positive(index=1, previous_phase_index=0, next_phase_index=0, dt=dt,
                                            time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[1],
                                            removal_at_phase_exit=removal_at_phase_exits[1],
                                            fixed_duration=fixed_durations[1], phase_duration=phase_durations[1],
                                            entry_function=entry_functions[1],
                                            entry_function_args=entry_functions_args[1],
                                            exit_function=exit_functions[1],
                                            exit_function_args=exit_functions_args[1],
                                            arrest_function=arrest_functions[1],
                                            arrest_function_args=arrest_functions_args[1],
                                            transition_to_next_phase=transitions_to_next_phase[1],
                                            transition_to_next_phase_args=transitions_to_next_phase_args[1],
                                            target_volume=target_volumes[1], volume=volumes[1],
                                            update_volume=update_volumes[1], update_volume_args=update_volumes_args[1],
                                            update_volume_rate=update_volume_rates[1],
                                            simulated_cell_volume=simulated_cell_volume)

        Ki67_negative = Phases.Ki67Negative(index=0, previous_phase_index=1, next_phase_index=1, dt=dt,
                                            time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[0],
                                            removal_at_phase_exit=removal_at_phase_exits[0],
                                            fixed_duration=fixed_durations[0], phase_duration=phase_durations[0],
                                            entry_function=entry_functions[0],
                                            entry_function_args=entry_functions_args[0],
                                            exit_function=exit_functions[0],
                                            exit_function_args=exit_functions_args[0],
                                            arrest_function=arrest_functions[0],
                                            arrest_function_args=arrest_functions_args[0],
                                            transition_to_next_phase=transitions_to_next_phase[0],
                                            transition_to_next_phase_args=transitions_to_next_phase_args[0],
                                            target_volume=target_volumes[0], volume=volumes[0],
                                            update_volume=update_volumes[0], update_volume_args=update_volumes_args[0],
                                            update_volume_rate=update_volume_rates[0],
                                            simulated_cell_volume=simulated_cell_volume)

        phases = [Ki67_negative, Ki67_positive]

        super().__init__(name=name, dt=dt, phases=phases, quiescent_phase=quiescent_phase, time_unit=time_unit)

    @staticmethod
    def _check_arguments(number_phases, cycle_name, target_volumes, division_at_phase_exits,
                         removal_at_phase_exits, fixed_durations, phase_durations, entry_functions,
                         entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
                         arrest_functions_args, transitions_to_next_phase, transitions_to_next_phase_args,
                         update_volumes, update_volumes_args, update_volume_rates):

        if len(target_volumes) != number_phases:
            raise ValueError(f"{cycle_name} has {number_phases} phases, {len(target_volumes)} target volumes defined")
        elif type(target_volumes) != list and type(target_volumes) != tuple:
            raise TypeError(f"`target_volumes` must be a list or tuple, got {type(target_volumes)}")

        if len(division_at_phase_exits) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(division_at_phase_exits)} division flags defined")
        elif type(division_at_phase_exits) != list and type(division_at_phase_exits) != tuple:
            raise TypeError(f"`division_at_phase_exits` must be a list or tuple, got {type(division_at_phase_exits)}")

        if len(removal_at_phase_exits) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(removal_at_phase_exits)} removal flags defined")
        elif type(removal_at_phase_exits) != list and type(removal_at_phase_exits) != tuple:
            raise TypeError(f"`removal_at_phase_exits` must be a list or tuple, got {type(removal_at_phase_exits)}")

        if type(fixed_durations) != list and type(fixed_durations) != tuple:
            raise TypeError(f"`fixed_durations` must be a list or tuple, got {type(fixed_durations)}")
        elif len(fixed_durations) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(fixed_durations)} fixed duration flags defined")

        if len(phase_durations) != number_phases:
            raise ValueError(f"{cycle_name} has {number_phases} phases, {len(phase_durations)} durations defined")
        elif type(phase_durations) != list and type(phase_durations) != tuple:
            raise TypeError(f"`phase_durations` must be a list or tuple, got {type(phase_durations)}")

        if len(entry_functions) != number_phases:
            raise ValueError(f"{cycle_name} has {number_phases} phases, {len(entry_functions)} entry functions defined")
        elif type(entry_functions) != list and type(entry_functions) != tuple:
            raise TypeError(f"`entry_functions` must be a list or tuple, got {type(entry_functions)}")

        if len(entry_functions_args) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(entry_functions_args)} entry functions args defined")
        elif type(entry_functions_args) != list and type(entry_functions_args) != tuple:
            raise TypeError(f"`entry_functions_args` must be a list or tuple, got {type(entry_functions_args)}")
        #
        if len(exit_functions) != number_phases:
            raise ValueError(f"{cycle_name} has {number_phases} phases, {len(exit_functions)} exit functions defined")
        elif type(exit_functions) != list and type(exit_functions) != tuple:
            raise TypeError(f"`exit_functions` must be a list or tuple, got {type(exit_functions)}")

        if len(exit_functions_args) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(exit_functions_args)} entry functions args defined")
        elif type(exit_functions_args) != list and type(exit_functions_args) != tuple:
            raise TypeError(f"`entry_functions_args` must be a list or tuple, got {type(exit_functions_args)}")
        #
        if len(arrest_functions) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(arrest_functions)} arrest functions defined")
        elif type(arrest_functions) != list and type(arrest_functions) != tuple:
            raise TypeError(f"`arrest_functions` must be a list or tuple, got {type(exit_functions)}")

        if len(arrest_functions_args) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(arrest_functions_args)} arrest functions args defined")
        elif type(arrest_functions_args) != list and type(arrest_functions_args) != tuple:
            raise TypeError(f"`arrest_functions_args` must be a list or tuple, got {type(arrest_functions_args)}")
        #
        if len(transitions_to_next_phase) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(transitions_to_next_phase)} transition functions defined")
        elif type(transitions_to_next_phase) != list and type(transitions_to_next_phase) != tuple:
            raise TypeError(
                f"`transitions_to_next_phase` must be a list or tuple, got {type(transitions_to_next_phase)}")

        if len(transitions_to_next_phase_args) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(transitions_to_next_phase_args)} transition functions args "
                f"defined")
        elif type(transitions_to_next_phase_args) != list and type(transitions_to_next_phase_args) != tuple:
            raise TypeError(
                f"`transitions_to_next_phase_args` must be a list or tuple, got {type(arrest_functions_args)}")

        #
        if len(update_volumes) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(update_volumes)} update volume functions defined")
        elif type(update_volumes) != list and type(update_volumes) != tuple:
            raise TypeError(f"`update_volumes` must be a list or tuple, got {type(update_volumes)}")

        if len(update_volumes_args) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(update_volumes_args)} update volume functions args "
                f"defined")
        elif type(update_volumes_args) != list and type(update_volumes_args) != tuple:
            raise TypeError(f"`update_volumes_args` must be a list or tuple, got {type(update_volumes_args)}")

        #
        if len(update_volume_rates) != number_phases:
            raise ValueError(
                f"{cycle_name} has {number_phases} phases, {len(update_volume_rates)} update volume rates defined")
        elif type(update_volume_rates) != list and type(update_volume_rates) != tuple:
            raise TypeError(f"`update_volume_rates` must be a list or tuple, got {type(update_volumes)}")


cycle_names = ["Simple Live", "Ki67 Basic"]


def get_cycle_by_name(name):
    if name not in cycle_names:
        raise ValueError(f"{name} is not a pre-defined cycle")

    if name == "Simple Live":
        return SimpleLiveCycle
    elif name == "Ki67 Basic":
        return Ki67Basic

    return Cycle


if __name__ == "__main__":
    print(cycle_names)

    test = Ki67Basic()
