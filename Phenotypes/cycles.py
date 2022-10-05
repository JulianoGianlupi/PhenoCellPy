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
        Time unit. TODO: Defines time convertions

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


class SimpleLiveCycle(Cycle):
    def __init__(self, time_unit: str = "min", name: str = "simple_live", dt=1):
        phases = [Phases.Phase(index=0, previous_phase_index=0, next_phase_index=0, dt=dt, time_unit=time_unit,
                               name="alive", division_at_phase_exit=True, phase_duration=60)]
        super().__init__(name=name, time_unit=time_unit, phases=phases, quiescent_phase=False, dt=dt)


class Ki67Basic(Cycle):
    def __init__(self, name="Ki67 Basic", dt=0.1, quiescent_phase=False):
        Ki67_positive = Phases.Ki67Negative(index=1, dt=dt, previous_phase_index=0, next_phase_index=0)
        Ki67_negative = Phases.Ki67Positive(index=0, dt=dt, previous_phase_index=1, next_phase_index=1)

        phases = [Ki67_negative, Ki67_positive]

        super().__init__(name=name, dt=dt, phases=phases, quiescent_phase=quiescent_phase)




