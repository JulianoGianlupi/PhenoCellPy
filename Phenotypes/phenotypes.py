import Phenotypes.phases as Phases
from numpy.random import randint


# todo:
#  - implement physicell's phenotypes
#  - documentation
#  - add biomodels ontology anotation
#  - switch how defalts are handled to be like phases.NecrosisLysed
#  - add copyright
#  - fix _check_arguments
#  - add cell heterogeneity. have (yet another) argument to set randomization to True/False (and possibly which distri-
#    tion to use, log-normal, normal, etc). If true, rates, durations, and volumes are slightly shuffled. Also have
#    an argument to set the variance


def _check_arguments(number_phases, cycle_name, division_at_phase_exits, removal_at_phase_exits, fixed_durations,
                     phase_durations, entry_functions, entry_functions_args, exit_functions, exit_functions_args,
                     arrest_functions, arrest_functions_args, transitions_to_next_phase, transitions_to_next_phase_args,
                     nuclear_biomass_change_rate, calcification_rate, cytoplasm_volume, cytoplasm_target_volume,
                     cytoplasm_target_fluid_fraction, nuclear_volume, nuclear_target_volume,
                     nuclear_target_fluid_fraction, calcified_fraction, cytoplasm_biomass_change_rate):
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
    if len(nuclear_biomass_change_rate) != number_phases:
        raise ValueError(
            f"{cycle_name} has {number_phases} phases, {len(nuclear_biomass_change_rate)} nuclear biomass change rates "
            f"defined")
    elif type(nuclear_biomass_change_rate) != list and type(nuclear_biomass_change_rate) != tuple:
        raise TypeError(
            f"`nuclear_biomass_change_rate` must be a list or tuple, got {type(nuclear_biomass_change_rate)}")

    #
    if len(calcification_rate) != number_phases:
        raise ValueError(
            f"{cycle_name} has {number_phases} phases, {len(calcification_rate)} calcification rates "
            f"defined")
    elif type(calcification_rate) != list and type(calcification_rate) != tuple:
        raise TypeError(
            f"`calcification_rate` must be a list or tuple, got {type(calcification_rate)}")

    #
    if len(cytoplasm_volume) != number_phases:
        raise ValueError(
            f"{cycle_name} has {number_phases} phases, {len(cytoplasm_volume)} cytoplasm volumes defined")
    elif type(cytoplasm_volume) != list and type(cytoplasm_volume) != tuple:
        raise TypeError(
            f"`cytoplasm_volume` must be a list or tuple, got {type(cytoplasm_volume)}")

    #
    if len(cytoplasm_target_volume) != number_phases:
        raise ValueError(
            f"{cycle_name} has {number_phases} phases, {len(cytoplasm_target_volume)} cytoplasm target volumes defined")
    elif type(cytoplasm_target_volume) != list and type(cytoplasm_target_volume) != tuple:
        raise TypeError(
            f"`cytoplasm_target_volume` must be a list or tuple, got {type(cytoplasm_target_volume)}")

        #
    if len(cytoplasm_target_fluid_fraction) != number_phases:
        raise ValueError(
            f"{cytoplasm_target_fluid_fraction} has {number_phases} phases, {len(cytoplasm_target_fluid_fraction)} "
            f"cytoplasm target fluid fractions defined")
    elif type(cytoplasm_target_fluid_fraction) != list and type(cytoplasm_target_fluid_fraction) != tuple:
        raise TypeError(
            f"`cytoplasm_target_fluid_fraction` must be a list or tuple, got {type(cytoplasm_target_fluid_fraction)}")

        #
    if len(nuclear_volume) != number_phases:
        raise ValueError(
            f"{cycle_name} has {number_phases} phases, {len(nuclear_volume)} nuclear volumes defined")
    elif type(nuclear_volume) != list and type(nuclear_volume) != tuple:
        raise TypeError(
            f"`nuclear_volume` must be a list or tuple, got {type(nuclear_volume)}")

        #
    if len(nuclear_target_volume) != number_phases:
        raise ValueError(
            f"{cycle_name} has {number_phases} phases, {len(nuclear_target_volume)} nuclear target volumes defined")
    elif type(nuclear_target_volume) != list and type(nuclear_target_volume) != tuple:
        raise TypeError(
            f"`nuclear_target_volume` must be a list or tuple, got {type(nuclear_target_volume)}")
    #
    if len(nuclear_target_fluid_fraction) != number_phases:
        raise ValueError(
            f"{cycle_name} has {number_phases} phases, {len(nuclear_target_fluid_fraction)} nuclear target fluid "
            f"fractions defined")
    elif type(nuclear_target_fluid_fraction) != list and type(nuclear_target_fluid_fraction) != tuple:
        raise TypeError(
            f"`nuclear_target_fluid_fraction` must be a list or tuple, got {type(nuclear_target_fluid_fraction)}")

    #
    if len(calcified_fraction) != number_phases:
        raise ValueError(
            f"{cycle_name} has {number_phases} phases, {len(calcified_fraction)} calcified fractions defined")
    elif type(calcified_fraction) != list and type(calcified_fraction) != tuple:
        raise TypeError(
            f"`calcified_fraction` must be a list or tuple, got {type(calcified_fraction)}")
        #
    if len(cytoplasm_biomass_change_rate) != number_phases:
        raise ValueError(
            f"{cycle_name} has {number_phases} phases, {len(cytoplasm_biomass_change_rate)} cytoplasm biomass change "
            f"rates defined")
    elif type(cytoplasm_biomass_change_rate) != list and type(cytoplasm_biomass_change_rate) != tuple:
        raise TypeError(
            f"`calcified_fraction` must be a list or tuple, got {type(cytoplasm_biomass_change_rate)}")


class Phenotype:
    """

    Base class to define a cell phenotype.

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

    starting_phase_index : int
        Index of which phase to start at. To start at a random phase set to -1

    quiescent_phase : :class:`Phases.Phase`, None, or False
        If false the cycle won't have a (stand-alone) quiescent phase defined. If None the default
        :class:`Phases.QuiescentPhase` will be used as the stand-alone quiescent phase. If it is a :class:`Phases.Phase`
        object it will be used as the stand-alone quiescent phase.

    current_phase : :class:`Phases.Phase`
        The current (active) phase of the cycle.

    time_in_phenotype : float
        Total time elapsed for the cycle

    """

    def __init__(self, name: str = "unnamed", dt: float = 1, time_unit: str = "min", phases: list = None,
                 quiescent_phase: Phases.Phase or False = None, starting_phase_index: int = 0):
        """

        :param name:
        :param dt:
        :param time_unit:
        :param phases:
        :param quiescent_phase:
        :param starting_phase_index:
        """
        # todo: add __init__ parameters for custom functions for each class
        # todo: add alias for self.current_phase.volume, i.e. property self.volume that fetches
        #  self.current_phase.volume
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
        if starting_phase_index is None:
            starting_phase_index = 0
        elif starting_phase_index == -1:  # random option
            starting_phase_index = randint(0, len(self.phases) + 1)

        self.current_phase = self.phases[starting_phase_index]
        self.time_in_phenotype = 0

    def time_step_phenotype(self):
        """
        Time-steps the cycle.

        Increments :attr:`time_in_cycle` by :attr:`dt`. Calls :func:`current_phase.time_step_phase`. If the phase time-
        step determines the cycle moves to the next phase (i.e., returns `True` for `next_phase`), calls
        :func:go_to_next_phase. If the phase time-step determines the cell exits the cell cycle and goes to quiescence
        (i.e., returns `True` for `quies`) calls :func:`go_to_quiescence`.

        :return: Flags (bool) for phase changing, cell death, and cell division
        :rtype: tuple of bool
        """
        # if we initialize a phenotype that has a initial phase with an entry function we should execute it immediately
        if not self.time_in_phenotype and self.current_phase.entry_function is not None:
            self.current_phase.entry_function(*self.current_phase.entry_function_args)

        self.time_in_phenotype += self.dt

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
            self.current_phase.entry_function(*self.current_phase.entry_function_args)
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

        # get the current cytoplasm, nuclear, calcified volumes
        cyto_solid = self.current_phase.volume.cytoplasm_solid
        cyto_fluid = self.current_phase.volume.cytoplasm_fluid

        nucl_solid = self.current_phase.volume.nuclear_solid
        nucl_fluid = self.current_phase.volume.nuclear_fluid

        calc_frac = self.current_phase.volume.calcified_fraction

        # get the target volumes

        target_cytoplasm_solid = self.current_phase.volume.cytoplasm_solid_target
        # target_cyto_fluid_frac = self.current_phase.volume.target_cytoplasm_fluid_fraction

        target_nuclear_solid = self.current_phase.volume.nuclear_solid_target
        # target_nucl_fluid_frac = self.current_phase.volume.target_nuclear_fluid_fraction

        target_fluid_fraction = self.current_phase.volume.target_fluid_fraction

        # set parameters of next phase

        self.phases[idx].volume.cytoplasm_solid = cyto_solid
        self.phases[idx].volume.cytoplasm_fluid = cyto_fluid

        self.phases[idx].volume.nuclear_solid = nucl_solid
        self.phases[idx].volume.nuclear_fluid = nucl_fluid

        self.phases[idx].volume.calcified_fraction = calc_frac

        self.phases[idx].volume.cytoplasm_solid_target = target_cytoplasm_solid

        self.phases[idx].volume.nuclear_solid_target = target_nuclear_solid

        self.phases[idx].volume.target_fluid_fraction = target_fluid_fraction

        # set phase

        self.current_phase = self.phases[idx]
        self.current_phase.time_in_phase = 0

        if self.current_phase.entry_function is not None:
            self.current_phase.entry_function(*self.current_phase.entry_function_args)

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

        # get the current cytoplasm, nuclear, calcified volumes
        cyto_solid = self.current_phase.volume.cytoplasm_solid
        cyto_fluid = self.current_phase.volume.cytoplasm_fluid

        nucl_solid = self.current_phase.volume.nuclear_solid
        nucl_fluid = self.current_phase.volume.nuclear_fluid

        calc_frac = self.current_phase.volume.calcified_fraction

        # get the target volumes

        # target_cytoplasm_solid = self.current_phase.volume.cytoplasm_solid_target
        # # target_cyto_fluid_frac = self.current_phase.volume.target_cytoplasm_fluid_fraction
        #
        # target_nuclear_solid = self.current_phase.volume.nuclear_solid_target
        # # target_nucl_fluid_frac = self.current_phase.volume.target_nuclear_fluid_fraction
        #
        # target_fluid_fraction = self.current_phase.volume.target_fluid_fraction

        # setting the quiescent phase volume parameters. As the cell is now quiescent it shouldn't want to change its
        # volume, so we set the targets to be the current measurements
        self.quiescent_phase.volume.cytoplasm_solid = cyto_solid
        self.quiescent_phase.volume.cytoplasm_fluid = cyto_fluid

        self.quiescent_phase.volume.nuclear_solid = nucl_solid
        self.quiescent_phase.volume.nuclear_fluid = nucl_fluid

        self.quiescent_phase.volume.nuclear_solid_target = nucl_solid
        self.quiescent_phase.volume.cytoplasm_solid_target = cyto_solid

        self.quiescent_phase.volume.calcified_fraction = calc_frac

        self.quiescent_phase.volume.target_fluid_fraction = (cyto_fluid + nucl_fluid) / (nucl_solid + nucl_fluid +
                                                                                         cyto_fluid + cyto_solid)

        # set the phase to be quiescent
        self.current_phase = self.quiescent_phase
        self.current_phase.time_in_phase = 0

    def __str__(self):
        phases = ""
        for p in self.phases:
            phases += f"{p}, "
        if len(phases) > 2:
            phases = phases[:-2]
        return f"{self.name} cycle, phases: {phases}"


class SimpleLiveCycle(Phenotype):
    """
        Simplest life cycle, it has only one phase. When "progressing" to the next phase it divides.
        """

    def __init__(self, time_unit: str = "min", name: str = "Simple Live", dt=1):
        phases = [
            Phases.Phase(index=0, previous_phase_index=0, next_phase_index=0, dt=dt, time_unit=time_unit, name="alive",
                         division_at_phase_exit=True, phase_duration=60 / 0.0432)]
        super().__init__(name=name, dt=dt, time_unit=time_unit, phases=phases, quiescent_phase=False)


class Ki67Basic(Phenotype):
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
                 simulated_cell_volume=None, cytoplasm_biomass_change_rate=(None, None),
                 nuclear_biomass_change_rate=(None, None), calcification_rate=(None, None), calcified_fraction=(0, 0),
                 target_fluid_fraction=(None, None), nuclear_fluid=(None, None), nuclear_solid=(None, None),
                 nuclear_solid_target=(None, None), cytoplasm_fluid=(None, None), cytoplasm_solid=(None, None),
                 cytoplasm_solid_target=(None, None), target_cytoplasm_to_nuclear_ratio=(None, None),
                 fluid_change_rate=(None, None)):
        # _check_arguments(2, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
        #                  entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
        #                  arrest_functions_args, transitions_to_next_phase, transitions_to_next_phase_args,
        #                  nuclear_biomass_change_rate, calcification_rate, cytoplasm_volume, cytoplasm_target_volume,
        #                  cytoplasm_target_fluid_fraction, nuclear_volume, nuclear_target_volume,
        #                  nuclear_target_fluid_fraction, calcified_fraction, cytoplasm_biomass_change_rate)

        Ki67_positive = Phases.Ki67Positive(index=1, previous_phase_index=0, next_phase_index=0, dt=dt,
                                            time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[1],
                                            removal_at_phase_exit=removal_at_phase_exits[1],
                                            fixed_duration=fixed_durations[1], phase_duration=phase_durations[1],
                                            entry_function=entry_functions[1],
                                            entry_function_args=entry_functions_args[1],
                                            exit_function=exit_functions[1], exit_function_args=exit_functions_args[1],
                                            arrest_function=arrest_functions[1],
                                            arrest_function_args=arrest_functions_args[1],
                                            transition_to_next_phase=transitions_to_next_phase[1],
                                            transition_to_next_phase_args=transitions_to_next_phase_args[1],
                                            simulated_cell_volume=simulated_cell_volume,
                                            cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[1],
                                            nuclear_biomass_change_rate=nuclear_biomass_change_rate[1],
                                            calcification_rate=calcification_rate[1],
                                            target_fluid_fraction=target_fluid_fraction[1],
                                            nuclear_fluid=nuclear_fluid[1], nuclear_solid=nuclear_solid[1],
                                            nuclear_solid_target=nuclear_solid_target[1],
                                            cytoplasm_fluid=cytoplasm_fluid[1], cytoplasm_solid=cytoplasm_solid[1],
                                            cytoplasm_solid_target=cytoplasm_solid_target[1],
                                            target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[1],
                                            calcified_fraction=calcified_fraction[1],
                                            fluid_change_rate=fluid_change_rate[1])

        Ki67_negative = Phases.Ki67Negative(index=0, previous_phase_index=1, next_phase_index=1, dt=dt,
                                            time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[0],
                                            removal_at_phase_exit=removal_at_phase_exits[0],
                                            fixed_duration=fixed_durations[0], phase_duration=phase_durations[0],
                                            entry_function=entry_functions[0],
                                            entry_function_args=entry_functions_args[0],
                                            exit_function=exit_functions[0], exit_function_args=exit_functions_args[0],
                                            arrest_function=arrest_functions[0],
                                            arrest_function_args=arrest_functions_args[0],
                                            transition_to_next_phase=transitions_to_next_phase[0],
                                            transition_to_next_phase_args=transitions_to_next_phase_args[0],
                                            simulated_cell_volume=simulated_cell_volume,
                                            cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[0],
                                            nuclear_biomass_change_rate=nuclear_biomass_change_rate[0],
                                            calcification_rate=calcification_rate[0],
                                            target_fluid_fraction=target_fluid_fraction[0],
                                            nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                                            nuclear_solid_target=nuclear_solid_target[0],
                                            cytoplasm_fluid=cytoplasm_fluid[0], cytoplasm_solid=cytoplasm_solid[0],
                                            cytoplasm_solid_target=cytoplasm_solid_target[0],
                                            target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                                            calcified_fraction=calcified_fraction[0],
                                            fluid_change_rate=fluid_change_rate[1])

        phases = [Ki67_negative, Ki67_positive]

        super().__init__(name=name, dt=dt, time_unit=time_unit, phases=phases, quiescent_phase=quiescent_phase)


class Ki67Advanced(Phenotype):

    def __init__(self, name="Ki67 Advanced", dt=0.1, time_unit="min", quiescent_phase=False,
                 division_at_phase_exits=(False, True, False), removal_at_phase_exits=(False, False, False),
                 fixed_durations=(False, True, True), phase_durations: list = (3.62 * 60, 13.0 * 60.0, 2.5 * 60),
                 entry_functions=(None, None, None), entry_functions_args=(None, None, None),
                 exit_functions=(None, None, None), exit_functions_args=(None, None, None),
                 arrest_functions=(None, None, None), arrest_functions_args=(None, None, None),
                 transitions_to_next_phase=(None, None, None),
                 transitions_to_next_phase_args: list = (None, None, None), simulated_cell_volume=None,
                 cytoplasm_biomass_change_rate=(None, None, None),
                 nuclear_biomass_change_rate=(None, None, None), calcification_rate=(None, None, None),
                 calcified_fraction=(0, 0, 0),
                 target_fluid_fraction=(None, None, None), nuclear_fluid=(None, None, None),
                 nuclear_solid=(None, None, None),
                 nuclear_solid_target=(None, None, None), cytoplasm_fluid=(None, None, None),
                 cytoplasm_solid=(None, None, None),
                 cytoplasm_solid_target=(None, None, None), target_cytoplasm_to_nuclear_ratio=(None, None, None),
                 fluid_change_rate=(None, None, None)):
        # _check_arguments(3, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
        #                  entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
        #                  arrest_functions_args, transitions_to_next_phase, transitions_to_next_phase_args,
        #                  nuclear_biomass_change_rate, calcification_rate, cytoplasm_volume, cytoplasm_target_volume,
        #                  cytoplasm_target_fluid_fraction, nuclear_volume, nuclear_target_volume,
        #                  nuclear_target_fluid_fraction, calcified_fraction, cytoplasm_biomass_change_rate)

        Ki67_negative = Phases.Ki67Negative(index=0, previous_phase_index=2, next_phase_index=1, dt=dt,
                                            time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[0],
                                            removal_at_phase_exit=removal_at_phase_exits[0],
                                            fixed_duration=fixed_durations[0], phase_duration=phase_durations[0],
                                            entry_function=entry_functions[0],
                                            entry_function_args=entry_functions_args[0],
                                            exit_function=exit_functions[0], exit_function_args=exit_functions_args[0],
                                            arrest_function=arrest_functions[0],
                                            arrest_function_args=arrest_functions_args[0],
                                            transition_to_next_phase=transitions_to_next_phase[0],
                                            transition_to_next_phase_args=transitions_to_next_phase_args[0],
                                            simulated_cell_volume=simulated_cell_volume,
                                            cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[0],
                                            nuclear_biomass_change_rate=nuclear_biomass_change_rate[0],
                                            calcification_rate=calcification_rate[0],
                                            target_fluid_fraction=target_fluid_fraction[0],
                                            nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                                            nuclear_solid_target=nuclear_solid_target[0],
                                            cytoplasm_fluid=cytoplasm_fluid[0], cytoplasm_solid=cytoplasm_solid[0],
                                            cytoplasm_solid_target=cytoplasm_solid_target[0],
                                            target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                                            calcified_fraction=calcified_fraction[0],
                                            fluid_change_rate=fluid_change_rate[0])

        Ki67_positive_pre = Phases.Ki67PositivePreMitotic(index=1, previous_phase_index=0, next_phase_index=2, dt=dt,
                                                          time_unit=time_unit,
                                                          division_at_phase_exit=division_at_phase_exits[1],
                                                          removal_at_phase_exit=removal_at_phase_exits[1],
                                                          fixed_duration=fixed_durations[1],
                                                          phase_duration=phase_durations[1],
                                                          entry_function=entry_functions[1],
                                                          entry_function_args=entry_functions_args[1],
                                                          exit_function=exit_functions[1],
                                                          exit_function_args=exit_functions_args[1],
                                                          arrest_function=arrest_functions[1],
                                                          arrest_function_args=arrest_functions_args[1],
                                                          transition_to_next_phase=transitions_to_next_phase[1],
                                                          transition_to_next_phase_args=transitions_to_next_phase_args[
                                                              1], simulated_cell_volume=simulated_cell_volume,
                                                          cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[
                                                              1],
                                                          nuclear_biomass_change_rate=nuclear_biomass_change_rate[1],
                                                          calcification_rate=calcification_rate[1],
                                                          target_fluid_fraction=target_fluid_fraction[1],
                                                          nuclear_fluid=nuclear_fluid[1],
                                                          nuclear_solid=nuclear_solid[1],
                                                          nuclear_solid_target=nuclear_solid_target[1],
                                                          cytoplasm_fluid=cytoplasm_fluid[1],
                                                          cytoplasm_solid=cytoplasm_solid[1],
                                                          cytoplasm_solid_target=cytoplasm_solid_target[1],
                                                          target_cytoplasm_to_nuclear_ratio=
                                                          target_cytoplasm_to_nuclear_ratio[1],
                                                          calcified_fraction=calcified_fraction[1],
                                                          fluid_change_rate=fluid_change_rate[1])

        Ki67_positive_post = Phases.Ki67PositivePostMitotic(index=2, previous_phase_index=1, next_phase_index=0, dt=dt,
                                                            time_unit=time_unit,
                                                            division_at_phase_exit=division_at_phase_exits[2],
                                                            removal_at_phase_exit=removal_at_phase_exits[2],
                                                            fixed_duration=fixed_durations[2],
                                                            phase_duration=phase_durations[2],
                                                            entry_function=entry_functions[2],
                                                            entry_function_args=entry_functions_args[2],
                                                            exit_function=exit_functions[2],
                                                            exit_function_args=exit_functions_args[2],
                                                            arrest_function=arrest_functions[2],
                                                            arrest_function_args=arrest_functions_args[2],
                                                            transition_to_next_phase=transitions_to_next_phase[2],
                                                            transition_to_next_phase_args=
                                                            transitions_to_next_phase_args[2],
                                                            simulated_cell_volume=simulated_cell_volume,
                                                            cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[
                                                                2],
                                                            nuclear_biomass_change_rate=nuclear_biomass_change_rate[2],
                                                            calcification_rate=calcification_rate[2],
                                                            target_fluid_fraction=target_fluid_fraction[2],
                                                            nuclear_fluid=nuclear_fluid[2],
                                                            nuclear_solid=nuclear_solid[2],
                                                            nuclear_solid_target=nuclear_solid_target[2],
                                                            cytoplasm_fluid=cytoplasm_fluid[2],
                                                            cytoplasm_solid=cytoplasm_solid[2],
                                                            cytoplasm_solid_target=cytoplasm_solid_target[2],
                                                            target_cytoplasm_to_nuclear_ratio=
                                                            target_cytoplasm_to_nuclear_ratio[2],
                                                            calcified_fraction=calcified_fraction[2],
                                                            fluid_change_rate=fluid_change_rate[2])
        phases = [Ki67_negative, Ki67_positive_pre, Ki67_positive_post]
        super().__init__(name=name, dt=dt, time_unit=time_unit, phases=phases, quiescent_phase=quiescent_phase)


class FlowCytometryBasic(Phenotype):
    """
    Basic flow cytometry model.

    Cell cycle consists of G0/G1 -> S -> G2/M -> G0/G1. Reference cycle lengths from
    https://www.ncbi.nlm.nih.gov/books/NBK9876/
    """

    def __init__(self, name="Flow Cytometry Basic", dt=0.1, time_unit="min", quiescent_phase=False,
                 division_at_phase_exits=(False, False, True), removal_at_phase_exits=(False, False, False),
                 fixed_durations=(False, False, False), phase_durations: list = (5.15 * 60, 8 * 60.0, 5 * 60),
                 entry_functions=(None, None, None), entry_functions_args=(None, None, None),
                 exit_functions=(None, None, None), exit_functions_args=(None, None, None),
                 arrest_functions=(None, None, None), arrest_functions_args=(None, None, None),
                 transitions_to_next_phase=(None, None, None),
                 transitions_to_next_phase_args: list = (None, None, None), simulated_cell_volume=None,
                 cytoplasm_biomass_change_rate=(None, None, None),
                 nuclear_biomass_change_rate=(None, None, None), calcification_rate=(None, None, None),
                 calcified_fraction=(0, 0, 0),
                 target_fluid_fraction=(None, None, None), nuclear_fluid=(None, None, None),
                 nuclear_solid=(None, None, None),
                 nuclear_solid_target=(None, None, None), cytoplasm_fluid=(None, None, None),
                 cytoplasm_solid=(None, None, None),
                 cytoplasm_solid_target=(None, None, None), target_cytoplasm_to_nuclear_ratio=(None, None, None),
                 fluid_change_rate=(None, None, None)):
        # _check_arguments(3, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
        #                  entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
        #                  arrest_functions_args, transitions_to_next_phase, transitions_to_next_phase_args,
        #                  nuclear_biomass_change_rate, calcification_rate, cytoplasm_volume, cytoplasm_target_volume,
        #                  cytoplasm_target_fluid_fraction, nuclear_volume, nuclear_target_volume,
        #                  nuclear_target_fluid_fraction, calcified_fraction, cytoplasm_biomass_change_rate)

        G0G1 = Phases.G0G1(dt=dt, time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[0],
                           removal_at_phase_exit=removal_at_phase_exits[0], fixed_duration=fixed_durations[0],
                           phase_duration=phase_durations[0], entry_function=entry_functions[0],
                           entry_function_args=entry_functions_args[0], exit_function=exit_functions[0],
                           exit_function_args=exit_functions_args[0], arrest_function=arrest_functions[0],
                           arrest_function_args=arrest_functions_args[0],
                           transition_to_next_phase=transitions_to_next_phase[0],
                           transition_to_next_phase_args=transitions_to_next_phase_args[0],
                           simulated_cell_volume=simulated_cell_volume,
                           cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[0],
                           nuclear_biomass_change_rate=nuclear_biomass_change_rate[0],
                           calcification_rate=calcification_rate[0], target_fluid_fraction=target_fluid_fraction[0],
                           nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                           nuclear_solid_target=nuclear_solid_target[0], cytoplasm_fluid=cytoplasm_fluid[0],
                           cytoplasm_solid=cytoplasm_solid[0], cytoplasm_solid_target=cytoplasm_solid_target[0],
                           target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                           calcified_fraction=calcified_fraction[0], fluid_change_rate=fluid_change_rate[0])

        S = Phases.S(dt=dt, time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[1],
                     removal_at_phase_exit=removal_at_phase_exits[1], fixed_duration=fixed_durations[1],
                     phase_duration=phase_durations[1], entry_function=entry_functions[1],
                     entry_function_args=entry_functions_args[1], exit_function=exit_functions[1],
                     exit_function_args=exit_functions_args[1], arrest_function=arrest_functions[1],
                     arrest_function_args=arrest_functions_args[1],
                     transition_to_next_phase=transitions_to_next_phase[1],
                     transition_to_next_phase_args=transitions_to_next_phase_args[1],
                     simulated_cell_volume=simulated_cell_volume,
                     cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[1],
                     nuclear_biomass_change_rate=nuclear_biomass_change_rate[1],
                     calcification_rate=calcification_rate[1], target_fluid_fraction=target_fluid_fraction[1],
                     nuclear_fluid=nuclear_fluid[1], nuclear_solid=nuclear_solid[1],
                     nuclear_solid_target=nuclear_solid_target[1], cytoplasm_fluid=cytoplasm_fluid[1],
                     cytoplasm_solid=cytoplasm_solid[1], cytoplasm_solid_target=cytoplasm_solid_target[1],
                     target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[1],
                     calcified_fraction=calcified_fraction[1], fluid_change_rate=fluid_change_rate[1])

        G2M = Phases.G2M(dt=dt, time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[2],
                         removal_at_phase_exit=removal_at_phase_exits[2], fixed_duration=fixed_durations[2],
                         phase_duration=phase_durations[2], entry_function=entry_functions[2],
                         entry_function_args=entry_functions_args[2], exit_function=exit_functions[2],
                         exit_function_args=exit_functions_args[2], arrest_function=arrest_functions[2],
                         arrest_function_args=arrest_functions_args[2],
                         transition_to_next_phase=transitions_to_next_phase[2],
                         transition_to_next_phase_args=transitions_to_next_phase_args[2],
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[2],
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate[2],
                         calcification_rate=calcification_rate[2], target_fluid_fraction=target_fluid_fraction[2],
                         nuclear_fluid=nuclear_fluid[2], nuclear_solid=nuclear_solid[2],
                         nuclear_solid_target=nuclear_solid_target[2], cytoplasm_fluid=cytoplasm_fluid[2],
                         cytoplasm_solid=cytoplasm_solid[2], cytoplasm_solid_target=cytoplasm_solid_target[2],
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[2],
                         calcified_fraction=calcified_fraction[2], fluid_change_rate=fluid_change_rate[2])

        phases = [G0G1, S, G2M]

        super().__init__(name=name, dt=dt, time_unit=time_unit, phases=phases, quiescent_phase=quiescent_phase)


class FlowCytometryAdvanced(Phenotype):
    """
    Basic flow cytometry model.

    Cell cycle consists of G0/G1 -> S -> G2 -> M -> G0/G1. Reference cycle lengths from
    https://www.ncbi.nlm.nih.gov/books/NBK9876/
    """

    def __init__(self, name="Flow Cytometry Advanced", dt=0.1, time_unit="min", quiescent_phase=False,
                 division_at_phase_exits=(False, False, False, True),
                 removal_at_phase_exits=(False, False, False, False), fixed_durations=(False, False, False, False),
                 phase_durations: list = (4.98 * 60, 8 * 60.0, 4 * 60, 1 * 60),
                 entry_functions=(None, None, None, None), entry_functions_args=(None, None, None, None),
                 exit_functions=(None, None, None, None), exit_functions_args=(None, None, None),
                 arrest_functions=(None, None, None, None), arrest_functions_args=(None, None, None, None),
                 transitions_to_next_phase=(None, None, None, None),
                 transitions_to_next_phase_args: list = (None, None, None, None), simulated_cell_volume=None,
                 cytoplasm_biomass_change_rate=(None, None, None, None),
                 nuclear_biomass_change_rate=(None, None, None, None), calcification_rate=(None, None, None, None),
                 calcified_fraction=(0, 0, 0, 0),
                 target_fluid_fraction=(None, None, None, None), nuclear_fluid=(None, None, None, None),
                 nuclear_solid=(None, None, None, None),
                 nuclear_solid_target=(None, None, None, None), cytoplasm_fluid=(None, None, None, None),
                 cytoplasm_solid=(None, None, None, None),
                 cytoplasm_solid_target=(None, None, None, None),
                 target_cytoplasm_to_nuclear_ratio=(None, None, None, None),
                 fluid_change_rate=(None, None, None, None)):
        # _check_arguments(4, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
        #                  entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
        #                  arrest_functions_args, transitions_to_next_phase, transitions_to_next_phase_args,
        #                  nuclear_biomass_change_rate, calcification_rate, cytoplasm_volume, cytoplasm_target_volume,
        #                  cytoplasm_target_fluid_fraction, nuclear_volume, nuclear_target_volume,
        #                  nuclear_target_fluid_fraction, calcified_fraction, cytoplasm_biomass_change_rate)

        G0G1 = Phases.G0G1(dt=dt, time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[0],
                           removal_at_phase_exit=removal_at_phase_exits[0], fixed_duration=fixed_durations[0],
                           phase_duration=phase_durations[0], entry_function=entry_functions[0],
                           entry_function_args=entry_functions_args[0], exit_function=exit_functions[0],
                           exit_function_args=exit_functions_args[0], arrest_function=arrest_functions[0],
                           arrest_function_args=arrest_functions_args[0],
                           transition_to_next_phase=transitions_to_next_phase[0],
                           transition_to_next_phase_args=transitions_to_next_phase_args[0],
                           simulated_cell_volume=simulated_cell_volume,
                           cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[0],
                           nuclear_biomass_change_rate=nuclear_biomass_change_rate[0],
                           calcification_rate=calcification_rate[0], target_fluid_fraction=target_fluid_fraction[0],
                           nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                           nuclear_solid_target=nuclear_solid_target[0], cytoplasm_fluid=cytoplasm_fluid[0],
                           cytoplasm_solid=cytoplasm_solid[0], cytoplasm_solid_target=cytoplasm_solid_target[0],
                           target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                           calcified_fraction=calcified_fraction[0], fluid_change_rate=fluid_change_rate[0])

        S = Phases.S(dt=dt, time_unit=time_unit, division_at_phase_exit=division_at_phase_exits[1],
                     removal_at_phase_exit=removal_at_phase_exits[1], fixed_duration=fixed_durations[1],
                     phase_duration=phase_durations[1], entry_function=entry_functions[1],
                     entry_function_args=entry_functions_args[1], exit_function=exit_functions[1],
                     exit_function_args=exit_functions_args[1], arrest_function=arrest_functions[1],
                     arrest_function_args=arrest_functions_args[1],
                     transition_to_next_phase=transitions_to_next_phase[1],
                     transition_to_next_phase_args=transitions_to_next_phase_args[1],
                     simulated_cell_volume=simulated_cell_volume,
                     cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[1],
                     nuclear_biomass_change_rate=nuclear_biomass_change_rate[1],
                     calcification_rate=calcification_rate[1], target_fluid_fraction=target_fluid_fraction[1],
                     nuclear_fluid=nuclear_fluid[1], nuclear_solid=nuclear_solid[1],
                     nuclear_solid_target=nuclear_solid_target[1], cytoplasm_fluid=cytoplasm_fluid[1],
                     cytoplasm_solid=cytoplasm_solid[1], cytoplasm_solid_target=cytoplasm_solid_target[1],
                     target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[1],
                     fluid_change_rate=fluid_change_rate[1])

        G2 = Phases.G0G1(index=2, previous_phase_index=1, next_phase_index=3, dt=dt, time_unit=time_unit, name="G2",
                         division_at_phase_exit=division_at_phase_exits[2],
                         removal_at_phase_exit=removal_at_phase_exits[2], fixed_duration=fixed_durations[2],
                         phase_duration=phase_durations[2], entry_function=entry_functions[2],
                         entry_function_args=entry_functions_args[2], exit_function=exit_functions[2],
                         exit_function_args=exit_functions_args[2], arrest_function=arrest_functions[2],
                         arrest_function_args=arrest_functions_args[2],
                         transition_to_next_phase=transitions_to_next_phase[2],
                         transition_to_next_phase_args=transitions_to_next_phase_args[2],
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[2],
                         nuclear_biomass_change_rate=nuclear_biomass_change_rate[2],
                         calcification_rate=calcification_rate[2], target_fluid_fraction=target_fluid_fraction[2],
                         nuclear_fluid=nuclear_fluid[2], nuclear_solid=nuclear_solid[2],
                         nuclear_solid_target=nuclear_solid_target[2], cytoplasm_fluid=cytoplasm_fluid[2],
                         cytoplasm_solid=cytoplasm_solid[2], cytoplasm_solid_target=cytoplasm_solid_target[2],
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[2],
                         fluid_change_rate=fluid_change_rate[2])

        M = Phases.G2M(index=3, previous_phase_index=2, next_phase_index=0, dt=dt, time_unit=time_unit, name="M",
                       division_at_phase_exit=division_at_phase_exits[3],
                       removal_at_phase_exit=removal_at_phase_exits[3], fixed_duration=fixed_durations[3],
                       phase_duration=phase_durations[3], entry_function=entry_functions[3],
                       entry_function_args=entry_functions_args[3], exit_function=exit_functions[3],
                       exit_function_args=exit_functions_args[3], arrest_function=arrest_functions[3],
                       arrest_function_args=arrest_functions_args[3],
                       transition_to_next_phase=transitions_to_next_phase[3],
                       transition_to_next_phase_args=transitions_to_next_phase_args[3],
                       simulated_cell_volume=simulated_cell_volume,
                       cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[3],
                       nuclear_biomass_change_rate=nuclear_biomass_change_rate[3],
                       calcification_rate=calcification_rate[3], target_fluid_fraction=target_fluid_fraction[3],
                       nuclear_fluid=nuclear_fluid[3], nuclear_solid=nuclear_solid[3],
                       nuclear_solid_target=nuclear_solid_target[3], cytoplasm_fluid=cytoplasm_fluid[3],
                       cytoplasm_solid=cytoplasm_solid[3], cytoplasm_solid_target=cytoplasm_solid_target[3],
                       target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[3],
                       fluid_change_rate=fluid_change_rate[3])

        phases = [G0G1, S, G2, M]

        super().__init__(name=name, dt=dt, time_unit=time_unit, phases=phases, quiescent_phase=quiescent_phase)


class ApoptosisStandard(Phenotype):
    """
    The standard apoptotic model
    """

    def __init__(self, name="Standard apoptosis model", dt=0.1, time_unit="min", quiescent_phase=False,
                 division_at_phase_exits=(False,), removal_at_phase_exits=(True,), fixed_durations=(True,),
                 phase_durations=(8.6 * 60,), entry_functions=(None,), entry_functions_args=(None,),
                 exit_functions=(None,), exit_functions_args=(None,), arrest_functions=(None,),
                 arrest_functions_args=(None,), transitions_to_next_phase=(None,),
                 transitions_to_next_phase_args=(None,), simulated_cell_volume=None,
                 cytoplasm_biomass_change_rate=(1 / 60,),
                 nuclear_biomass_change_rate=(0.35 / 60,), calcification_rate=(0,),
                 calcified_fraction=(0,),
                 target_fluid_fraction=(None,), nuclear_fluid=(None,), nuclear_solid=(None,),
                 nuclear_solid_target=(None,), cytoplasm_fluid=(None,), cytoplasm_solid=(None,),
                 cytoplasm_solid_target=(None,), target_cytoplasm_to_nuclear_ratio=(None,),
                 fluid_change_rate=(None,)):
        # _check_arguments(1, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
        #                  entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
        #                  arrest_functions_args, transitions_to_next_phase, transitions_to_next_phase_args,
        #                  nuclear_biomass_change_rate, calcification_rate, cytoplasm_volume, cytoplasm_target_volume,
        #                  cytoplasm_target_fluid_fraction, nuclear_volume, nuclear_target_volume,
        #                  nuclear_target_fluid_fraction, calcified_fraction, cytoplasm_biomass_change_rate)

        apopto = Phases.Apoptosis(index=0, previous_phase_index=0, next_phase_index=0, dt=dt, time_unit=time_unit,
                                  name="Apoptosis", division_at_phase_exit=division_at_phase_exits[0],
                                  removal_at_phase_exit=removal_at_phase_exits[0], fixed_duration=fixed_durations[0],
                                  phase_duration=phase_durations[0], entry_function=entry_functions[0],
                                  entry_function_args=entry_functions_args[0], exit_function=exit_functions[0],
                                  exit_function_args=exit_functions_args[0], arrest_function=arrest_functions[0],
                                  arrest_function_args=arrest_functions_args[0],
                                  transition_to_next_phase=transitions_to_next_phase[0],
                                  transition_to_next_phase_args=transitions_to_next_phase_args[0],
                                  simulated_cell_volume=simulated_cell_volume,
                                  cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[0],
                                  nuclear_biomass_change_rate=nuclear_biomass_change_rate[0],
                                  calcification_rate=calcification_rate[0],
                                  target_fluid_fraction=target_fluid_fraction[0], nuclear_fluid=nuclear_fluid[0],
                                  nuclear_solid=nuclear_solid[0], nuclear_solid_target=nuclear_solid_target[0],
                                  cytoplasm_fluid=cytoplasm_fluid[0], cytoplasm_solid=cytoplasm_solid[0],
                                  cytoplasm_solid_target=cytoplasm_solid_target[0],
                                  target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                                  calcified_fraction=calcified_fraction[0],
                                  fluid_change_rate=fluid_change_rate[0])

        # a phase to help lyse the simulated cell, shouldn't do anything
        # debris = Phases.Phase(index=1, previous_phase_index=0, next_phase_index=1, dt=dt, time_unit=time_unit,
        #                       name="Debris", division_at_phase_exit=False, removal_at_phase_exit=True,
        #                       fixed_duration=True, phase_duration=1e6)

        phases = [apopto]

        super().__init__(name=name, dt=dt, time_unit=time_unit, phases=phases, quiescent_phase=quiescent_phase)


class NecrosisStandard(Phenotype):
    """
    Standard Necrosis model
    """

    def __init__(self, name="Standard necrosis model", dt=0.1, time_unit="min", quiescent_phase=False,
                 division_at_phase_exits=(False, False), removal_at_phase_exits=(False, True),
                 fixed_durations=(False, True),
                 phase_durations=(None, None), entry_functions=(None, None), entry_functions_args=(None, None),
                 exit_functions=(None, None), exit_functions_args=(None, None), arrest_functions=(None, None),
                 arrest_functions_args=(None, None), transitions_to_next_phase=(None, None),
                 transitions_to_next_phase_args=(None, None), simulated_cell_volume=None,
                 cytoplasm_biomass_change_rate=(None, None),
                 nuclear_biomass_change_rate=(None, None),
                 calcification_rate=(None, None),
                 calcified_fraction=(0, 0),
                 target_fluid_fraction=(None, None), nuclear_fluid=(None, None), nuclear_solid=(None, None),
                 nuclear_solid_target=(None, None), cytoplasm_fluid=(None, None), cytoplasm_solid=(None, None),
                 cytoplasm_solid_target=(None, None), target_cytoplasm_to_nuclear_ratio=(None, None),
                 fluid_change_rate=(None, None)):

        necro_swell = Phases.NecrosisSwell(index=0, previous_phase_index=0, next_phase_index=1, dt=dt,
                                           time_unit=time_unit,
                                           division_at_phase_exit=division_at_phase_exits[0],
                                           removal_at_phase_exit=removal_at_phase_exits[0],
                                           fixed_duration=fixed_durations[0],
                                           phase_duration=phase_durations[0], entry_function=entry_functions[0],
                                           entry_function_args=entry_functions_args[0], exit_function=exit_functions[0],
                                           exit_function_args=exit_functions_args[0],
                                           arrest_function=arrest_functions[0],
                                           arrest_function_args=arrest_functions_args[0],
                                           transition_to_next_phase=transitions_to_next_phase[0],
                                           transition_to_next_phase_args=transitions_to_next_phase_args[0],
                                           simulated_cell_volume=simulated_cell_volume,
                                           cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[0],
                                           nuclear_biomass_change_rate=nuclear_biomass_change_rate[0],
                                           calcification_rate=calcification_rate[0],
                                           target_fluid_fraction=target_fluid_fraction[0],
                                           nuclear_fluid=nuclear_fluid[0],
                                           nuclear_solid=nuclear_solid[0], nuclear_solid_target=nuclear_solid_target[0],
                                           cytoplasm_fluid=cytoplasm_fluid[0], cytoplasm_solid=cytoplasm_solid[0],
                                           cytoplasm_solid_target=cytoplasm_solid_target[0],
                                           target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                                           calcified_fraction=calcified_fraction[0],
                                           fluid_change_rate=fluid_change_rate[0])

        necro_lysed = Phases.NecrosisLysed(index=1, previous_phase_index=0, next_phase_index=1, dt=dt,
                                           time_unit=time_unit,
                                           division_at_phase_exit=division_at_phase_exits[1],
                                           removal_at_phase_exit=removal_at_phase_exits[1],
                                           fixed_duration=fixed_durations[1],
                                           phase_duration=phase_durations[1], entry_function=entry_functions[1],
                                           entry_function_args=entry_functions_args[1], exit_function=exit_functions[1],
                                           exit_function_args=exit_functions_args[1],
                                           arrest_function=arrest_functions[1],
                                           arrest_function_args=arrest_functions_args[1],
                                           transition_to_next_phase=transitions_to_next_phase[1],
                                           transition_to_next_phase_args=transitions_to_next_phase_args[1],
                                           simulated_cell_volume=simulated_cell_volume,
                                           cytoplasm_biomass_change_rate=cytoplasm_biomass_change_rate[1],
                                           nuclear_biomass_change_rate=nuclear_biomass_change_rate[1],
                                           calcification_rate=calcification_rate[1],
                                           target_fluid_fraction=target_fluid_fraction[1],
                                           nuclear_fluid=nuclear_fluid[1],
                                           nuclear_solid=nuclear_solid[1], nuclear_solid_target=nuclear_solid_target[1],
                                           cytoplasm_fluid=cytoplasm_fluid[1], cytoplasm_solid=cytoplasm_solid[1],
                                           cytoplasm_solid_target=cytoplasm_solid_target[1],
                                           target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[1],
                                           calcified_fraction=calcified_fraction[1],
                                           fluid_change_rate=fluid_change_rate[1])

        phases = [necro_swell, necro_lysed]

        super().__init__(name=name, dt=dt, time_unit=time_unit, phases=phases, quiescent_phase=quiescent_phase)

        return


cycle_names = ["Simple Live", "Ki67 Basic"]


def get_phenotype_by_name(name):  # todo: keep updated
    if name not in cycle_names:
        raise ValueError(f"{name} is not a pre-defined cycle")

    if name == "Simple Live":
        return SimpleLiveCycle
    elif name == "Ki67 Basic":
        return Ki67Basic


    return Phenotype


if __name__ == "__main__":
    print(cycle_names)

    test = Ki67Basic()
