from numpy import exp
from numpy.random import uniform

from Phenotypes.cell_volume import CellVolumes


# todo: change args handling to also accept tuples

class Phase:
    """
    Base class to define phases of a cell phenotype.

    This base class is inherited by all phases. It defines some methods to time-step the Phase model, to transition to
    the next phase, and what should be done by the phenotype model when entering/exiting a phase. Initializes the cell
    volume model.

    Methods:
    -------

    time_step_phase()
        Time steps the phase. Returns a tuple (did the cell transition to the next phase, did the cell enter
        quiescence). See `time_step_phase`'s documentation for further explanation.

    transition_to_next_phase(*args)
        One of the default transition functions (`_transition_to_next_phase_deterministic`,
        `_transition_to_next_phase_stochastic`) or a user defined function. If user defined, it will be called with
        `transition_to_next_phase_args` as args. Must return a bool denoting if the transition occurs or not.

    _transition_to_next_phase_deterministic()
        Default deterministic transition function. Returns `time_in_phase > phase_duration`

    _transition_to_next_phase_stochastic()
        Default stochastic transition function. Probability of transition depends on `dt` and `phase_duration`

    entry_function(*args)
        Optional function to be executed upon entering this phase. Some pre-build Phases have their own entry function
        already defined. It gets called using attribute `entry_function_args`. Must have no return

    exit_function(*args)
        Optional function to be executed just before exiting this phase. Some pre-build Phases have their own exit
        function already defined. It gets called using attribute `exit_function_args`. Must have no return

    arrest_function(*args)
        Optional function that returns true if the cell should exit the cell cycle and enter quiescence

    update_volume()
        Function to update the volume of the cell. Calls the cell volume submodel (class:CellVolume) `update_volume`
        function

    _double_target_volume()
        Function that doubles all the target volumes of subclass class:CellVolume. Is used by some inheriting Phases as
        the entry function

    Parameters:
    -----------

    :param index: Index of this phase in the list of phases that forms a phenotype
    :type index: int

    :param previous_phase_index: Index of the phase preceding this phase in the list of phases that forms a
    phenotype
    :type previous_phase_index: int

    :param next_phase_index: Index of the phase proceding this phase in the list of phases that forms a phenotype
    :type next_phase_index: int

    :param dt: Time step duration in units of `time_unit`. `dt > 0`
    :type dt: float

    :param time_unit: What time units are used by the model (e.g., minutes, hours, days)
    :type time_unit: str

    :param name: Descriptive name of this phase (e.g., S, G, M, necrotic swelling)
    :type name: str

    :param division_at_phase_exit: Boolean to define if the simulated cell should divide when switching phases.
    :type division_at_phase_exit: bool

    :param removal_at_phase_exit: Boolean defining if the simulated cell should be removed (i.e., it dies, is
    killed, leaves the simulated domain) when leaving this phase
    :type removal_at_phase_exit: bool

    :param fixed_duration: Boolean seting the transition from this phase to the next to be deterministic (True) or
    stochastic (False)
    :type fixed_duration: bool

    :param phase_duration: Time duration, in units of `time_unit`, of this phase. In the case of the stochastic
    transition, the transition rate will be `dt/phase_duration`. `phase_duration > 0`
    :type phase_duration: float

    :param entry_function: Function that gets called immediately when entering this phase. Some phases have a
    default `entry_function` defined, such as Apoptosis. Must be an *args function
    :type entry_function: function

    :param entry_function_args: Args for `entry_function`
    :type entry_function_args: list or tuple

    :param exit_function: Function that gets called just before exiting this phase. Some phases have a
    default `exit_function` defined.
    :type exit_function: function

    :param exit_function_args: Args for `exit_function`
    :type exit_function_args: list or tuple

    :param arrest_function: Function that return if the cell should exit the current phase early and the cell cyle
    and enter quiescence
    :type arrest_function: function

    :param arrest_function_args: Args for `arrest_function`
    :type arrest_function_args: list or tuple

    :param transition_to_next_phase: Default or custom function that returns if the cell should advance to the next
    phase in the phenotype. If left as None the phase will pick either `_transition_to_next_phase_deterministic`
    or `_transition_to_next_phase_stochastic` depending on the value of `fixed_duration`.
    :type transition_to_next_phase: function

    :param transition_to_next_phase_args: Args for `transition_to_next_phase`
    :type transition_to_next_phase_args: list or tuple

    :param simulated_cell_volume: Volume of the simulated cell (e.g., a CompuCell3D or Tissue Forge cell)
    :type simulated_cell_volume: float

    :param cytoplasm_biomass_change_rate: Change rate for the cytoplasmic volume. volume/`time_unit` units.
    `cytoplasm_biomass_change_rate` >= 0. Passed to the `CellVolume` attribute class
    :type cytoplasm_biomass_change_rate: float

    :param nuclear_biomass_change_rate: Change rate for the nuclear volume. volume/`time_unit` units.
    `nuclear_biomass_change_rate` >= 0. Passed to the `CellVolume` attribute class
    :type nuclear_biomass_change_rate: float

    :param calcification_rate: Rate of calcification of the cell. volume/`time_unit` units. `calcification_rate` >= 0
    Passed to the `CellVolume` attribute class
    :type calcification_rate: float

    :param target_fluid_fraction: Fraction of the cell volume it will attempt to keep as fluid.
    0 <= `target_fluid_fraction` <= 1. Passed to the `CellVolume` attribute class
    :type target_fluid_fraction: float

    :param nuclear_fluid: Initial volume of the fluid part of the nucleus. Passed to the `CellVolume` attribute class
    :type nuclear_fluid: float

    :param nuclear_solid: Initial volume of the solid part of the nucleus. Passed to the `CellVolume` attribute class
    :type nuclear_solid: float

    :param nuclear_solid_target: Nuclear solid volume the volume model will tend towards. `nuclear_solid_target`>=0.
    Passed to the `CellVolume` attribute class
    :type nuclear_solid_target: float

    :param cytoplasm_fluid: Initial volume of the fluid part of the cytoplasm. Passed to the `CellVolume` attribute
    class
    :type cytoplasm_fluid: float

    :param cytoplasm_solid: Initial volume of the solid part of the cytoplasm. Passed to the `CellVolume` attribute
    class
    :type cytoplasm_solid: float

    :param cytoplasm_solid_target: Cytoplasm solid volume the volume model will tend towards.
    `cytoplasm_solid_target`>=0. Passed to the `CellVolume` attribute class
    :type cytoplasm_solid_target: float

    :param target_cytoplasm_to_nuclear_ratio: Nuclear to cytoplasmic volumes ratio the volume model will tend
    towards. `target_cytoplasm_to_nuclear_ratio`>0
    :type target_cytoplasm_to_nuclear_ratio: float

    :param calcified_fraction: Initial fraction of the cell that is calcified. 0<=`calcified_fraction`<=1
    :type calcified_fraction: float

    :param fluid_change_rate: Rate of change of the cell fluid part. volume/`time_unit` units.
    `fluid_change_rate` >= 0. Passed to the `CellVolume` attribute class
    :type fluid_change_rate: float

    :param relative_rupture_volume: Proportion of the initial volume that causes cell lysis
    :type relative_rupture_volume: float

    Attributes
    ----------

    time_in_phase : float
        Time spent in this phase

    volume : class:cell_volume.CellVolumes
        Cell volume submodel

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
                 calcified_fraction=None, fluid_change_rate=None, relative_rupture_volume=None):
        """

        :param index: Index of this phase in the list of phases that forms a phenotype
        :type index: int

        :param previous_phase_index: Index of the phase preceding this phase in the list of phases that forms a
        phenotype
        :type previous_phase_index: int

        :param next_phase_index: Index of the phase proceding this phase in the list of phases that forms a phenotype
        :type next_phase_index: int

        :param dt: Time step duration in units of `time_unit`. `dt > 0`
        :type dt: float

        :param time_unit: What time units are used by the model (e.g., minutes, hours, days)
        :type time_unit: str

        :param name: Descriptive name of this phase (e.g., S, G, M, necrotic swelling)
        :type name: str

        :param division_at_phase_exit: Boolean to define if the simulated cell should divide when switching phases.
        :type division_at_phase_exit: bool

        :param removal_at_phase_exit: Boolean defining if the simulated cell should be removed (i.e., it dies, is
        killed, leaves the simulated domain) when leaving this phase
        :type removal_at_phase_exit: bool

        :param fixed_duration: Boolean seting the transition from this phase to the next to be deterministic (True) or
        stochastic (False)
        :type fixed_duration: bool

        :param phase_duration: Time duration, in units of `time_unit`, of this phase. In the case of the stochastic
        transition, the transition rate will be `dt/phase_duration`. `phase_duration > 0`
        :type phase_duration: float

        :param entry_function: Function that gets called immediately when entering this phase. Some phases have a
        default `entry_function` defined, such as Apoptosis. Must be an *args function
        :type entry_function: function

        :param entry_function_args: Args for `entry_function`
        :type entry_function_args: list or tuple

        :param exit_function: Function that gets called just before exiting this phase. Some phases have a
        default `exit_function` defined.
        :type exit_function: function

        :param exit_function_args: Args for `exit_function`
        :type exit_function_args: list or tuple

        :param arrest_function: Function that return if the cell should exit the current phase early and the cell cyle
        and enter quiescence
        :type arrest_function: function

        :param arrest_function_args: Args for `arrest_function`
        :type arrest_function_args: list or tuple

        :param transition_to_next_phase: Default or custom function that returns if the cell should advance to the next
        phase in the phenotype. If left as None the phase will pick either `_transition_to_next_phase_deterministic`
        or `_transition_to_next_phase_stochastic` depending on the value of `fixed_duration`.
        :type transition_to_next_phase: function

        :param transition_to_next_phase_args: Args for `transition_to_next_phase`
        :type transition_to_next_phase_args: list or tuple

        :param simulated_cell_volume: Volume of the simulated cell (e.g., a CompuCell3D or Tissue Forge cell)
        :type simulated_cell_volume: float

        :param cytoplasm_biomass_change_rate: Change rate for the cytoplasmic volume. volume/`time_unit` units.
        `cytoplasm_biomass_change_rate` >= 0. Passed to the `CellVolume` attribute class
        :type cytoplasm_biomass_change_rate: float

        :param nuclear_biomass_change_rate: Change rate for the nuclear volume. volume/`time_unit` units.
        `nuclear_biomass_change_rate` >= 0. Passed to the `CellVolume` attribute class
        :type nuclear_biomass_change_rate: float

        :param calcification_rate: Rate of calcification of the cell. volume/`time_unit` units. `calcification_rate` >= 0
        Passed to the `CellVolume` attribute class
        :type calcification_rate: float

        :param target_fluid_fraction: Fraction of the cell volume it will attempt to keep as fluid.
        0 <= `target_fluid_fraction` <= 1. Passed to the `CellVolume` attribute class
        :type target_fluid_fraction: float

        :param nuclear_fluid: Initial volume of the fluid part of the nucleus. Passed to the `CellVolume` attribute class
        :type nuclear_fluid: float

        :param nuclear_solid: Initial volume of the solid part of the nucleus. Passed to the `CellVolume` attribute class
        :type nuclear_solid: float

        :param nuclear_solid_target: Nuclear solid volume the volume model will tend towards. `nuclear_solid_target`>=0.
        Passed to the `CellVolume` attribute class
        :type nuclear_solid_target: float

        :param cytoplasm_fluid: Initial volume of the fluid part of the cytoplasm. Passed to the `CellVolume` attribute
        class
        :type cytoplasm_fluid: float

        :param cytoplasm_solid: Initial volume of the solid part of the cytoplasm. Passed to the `CellVolume` attribute
        class
        :type cytoplasm_solid: float

        :param cytoplasm_solid_target: Cytoplasm solid volume the volume model will tend towards.
        `cytoplasm_solid_target`>=0. Passed to the `CellVolume` attribute class
        :type cytoplasm_solid_target: float

        :param target_cytoplasm_to_nuclear_ratio: Nuclear to cytoplasmic volumes ratio the volume model will tend
        towards. `target_cytoplasm_to_nuclear_ratio`>0
        :type target_cytoplasm_to_nuclear_ratio: float

        :param calcified_fraction: Initial fraction of the cell that is calcified. 0<=`calcified_fraction`<=1
        :type calcified_fraction: float

        :param fluid_change_rate: Rate of change of the cell fluid part. volume/`time_unit` units.
        `fluid_change_rate` >= 0. Passed to the `CellVolume` attribute class
        :type fluid_change_rate: float

        :param relative_rupture_volume: Proportion of the initial volume that causes cell lysis
        :type relative_rupture_volume: float
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
        if self.exit_function is not None and not (type(self.exit_function_args) == list or
                                                   type(self.exit_function_args) == tuple):
            raise TypeError("Exit function defined but no args given. Was expecting "
                            f"'exit_function_args' to be a list or tupple, got {type(exit_function_args)}.")

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
                                  calcified_fraction=calcified_fraction,
                                  relative_rupture_volume=relative_rupture_volume)

    def update_volume(self):
        """
        Calls the cell volume submodel :function:`CellVolumes.update_volume`. Passes the current phase volume change
        rates as well as the timestep to it.

        :return: No return
        """
        self.volume.update_volume(self.dt, self.fluid_change_rate, self.nuclear_biomass_change_rate,
                                  self.cytoplasm_biomass_change_rate, self.calcification_rate)

    def _transition_to_next_phase_stochastic(self, *none):
        """
        Default stochastic phase transition function.

        Calculates a Poisson probability based on dt and self.phase_duration (p=1-exp(-dt/phase_duration), rolls a
        random number, and returns random number < probability.

        :param none: Not used. Placeholder in case of user defined function with args
        :return: bool. random number < probability of transition
        """

        # the approximation 1-exp(-x) ~ x can be used. That approximation has a difference of 0.005 at x=0.1, which I'd
        # find acceptable. TODO: implement a check on self.dt / self.phase_duration, if it is < .1 use the approximation

        prob = 1 - exp(-self.dt / self.phase_duration)
        return uniform() < prob

    def _transition_to_next_phase_deterministic(self, *none):
        """
        Default deterministic phase transition function.

        If the time spent in this phase is greater than the phase duration, go to the next phase.

        :param none: Not used. Placeholder in case of user defined function with args
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
        """

        Doubles the cell volume submodel (:class:`Phenotypes.cell_volume`) target volumes. Used by several cell cycle
        models to double the cell volume before mitosis

        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return: No return
        """
        self.volume.nuclear_solid_target *= 2
        self.volume.cytoplasm_solid_target *= 2

    def _halve_target_volume(self, *none):
        """

        Halves the cell volume submodel (:class:`Phenotypes.cell_volume`) target volumes. Used by several cell cycle
        models to halve the cell volume after mitosis

        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return: No return
        """
        self.volume.cytoplasm_solid_target /= 2
        self.volume.nuclear_solid_target /= 2

    def __str__(self):
        return f"{self.name} phase"


class QuiescentPhase(Phase):
    """
    Default Quiescent Phase. Inherits :class:`Phase`

    This quiescent phase class is meant to be "outside" whatever phenotype progression is being used.

    """

    def __init__(self, index: int = 9999, previous_phase_index: int = None, next_phase_index: int = 0, dt: float = None,
                 time_unit: str = "min", name: str = "quiescent", division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 4.59 * 60,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate=0, nuclear_biomass_change_rate=0,
                 calcification_rate=0, target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None,
                 nuclear_solid_target=None, cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None,
                 relative_rupture_volume=None):
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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)
        return


class Ki67Negative(Phase):
    """
    Inherits :class:`Phase`. Defines Ki 67- quiescent phase.

    This is a quiescent phenotype for cells that are replicating. Ki67 is a protein marker associated with proliferation.
    Transition to the next phase is set to be stochastic (the phase does not use a fixed duration) by default. Default
    expected phase duration is 4.59h, the phase transition rate is, therefore, dt/4.59 1/h.
    This phase does not calcify the cell. The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939
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
                 calcified_fraction=None, fluid_change_rate=None, relative_rupture_volume=None):
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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)


class Ki67Positive(Phase):
    """

    Inherits :class:`Phase`. Defines Ki 67+ proliferating phase.

    This is a proliferating phenotype for cells that are replicating. Ki67 is a protein marker associated with prolife-
    ration. Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by default.
    Default phase duration is 15.5h. By default, if no user defined custom entry function is defined (i.e.,
    `entry_function=None`), this phase will set its entry function to be :class:`Phase._double_target_volume`.
    By default, will set the volume change rates to be [change in volume]/[phase duration]. This phase does not calcify
    the cell.

    The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939

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
                 calcified_fraction=None, fluid_change_rate=None, relative_rupture_volume=None):

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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)


class Ki67PositivePreMitotic(Ki67Positive):
    """

    Inherits :class:`Ki67Positive`. Defines Ki 67+ pre-mitotic proliferating phase. Only difference to
    :class:`Ki67Positive` is the phase length.

    This is a proliferating phenotype for cells that are replicating. Ki67 is a protein marker associated with prolife-
    ration. Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by default.
    Default phase duration is 13h. By default, if no user defined custom entry function is defined (i.e.,
    `entry_function=None`), this phase will set its entry function to be :class:`Phase._double_target_volume`

    The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939

    """

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
                 calcified_fraction=None, fluid_change_rate=None, relative_rupture_volume=None):
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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)


class Ki67PositivePostMitotic(Phase):
    """
    Inherits :class:`Phase`. Defines Ki 67+ post-mitotic phase, it represents the cell's reorganization.

    This is a rest phenotype for cells that are replicating. Ki67 is a protein marker associated with proliferation.
    Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by default. Default
    phase duration is 2.5h. By default, if no user defined custom entry function is defined (i.e.,
    `entry_function=None`), this phase will set its entry function to be
    :class:`Ki67PositivePostMitotic._standard_Ki67_positive_postmit_entry_function`, which calls
    :class:`Phase._halve_target_volume`.

    The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939
    """
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
                 calcified_fraction=None, fluid_change_rate=None, relative_rupture_volume=None):

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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)

    def _standard_Ki67_positive_postmit_entry_function(self, *args):
        """
        Calls :class:`Phase._halve_target_volume`.

        :param args: Not used
        :return:
        """
        self._halve_target_volume(*args)


class G0G1(Phase):
    """
    Inherits :class:`Phase`. Defines G0/G1 phase, it more representative of the quiescent phase than the first growth
    phase.

    This is a quiescent phenotype for cells that are replicating. Transition to the next phase is set to be stochastic
    (the phase does not use a fixed duration) by default. Default
    expected phase duration is 5.15h, the phase transition rate is, therefore, dt/5.15 1/h.
    This phase does not calcify the cell. Reference phase duration from https://www.ncbi.nlm.nih.gov/books/NBK9876/
    """
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
                 calcified_fraction=None, fluid_change_rate=None, relative_rupture_volume=None):
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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)


class S(Phase):
    """
    Inherits :class:`Phase`. Defines S phase, it more representative of the growth phase than the inter-growth rest.

    This is a growth phenotype for cells that are replicating. Transition to the next phase is set to be stochastic
    (the phase does not use a fixed duration) by default. Default
    expected phase duration is 8h, the phase transition rate is, therefore, dt/8 1/h. By default, will set the volume
    change rates to be [change in volume]/[phase duration]. This phase does not calcify the cell. Reference phase
    duration from https://www.ncbi.nlm.nih.gov/books/NBK9876/
    """
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
                 calcified_fraction=None, fluid_change_rate=None, relative_rupture_volume=None):

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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)


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
                 calcified_fraction=None, fluid_change_rate=None, relative_rupture_volume=None):
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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)


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

        if entry_function is None:
            entry_function = self._standard_apoptosis_entry
            entry_function_args = [None]

        if unlysed_fluid_change_rate is not None:
            self.unlysed_fluid_change_rate = unlysed_fluid_change_rate
        else:
            self.unlysed_fluid_change_rate = 3 / 60

        if lysed_fluid_change_rate is not None:
            self.lysed_fluid_change_rate = lysed_fluid_change_rate
        else:
            self.lysed_fluid_change_rate = 0

        if relative_rupture_volume is None:
            self.relative_rupture_volume = 2
        else:
            self.relative_rupture_volume = relative_rupture_volume

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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)

    def _standard_apoptosis_entry(self, *none):

        # shrink cell
        self.volume.target_fluid_fraction = 0
        self.volume.cytoplasm_solid_target = 0
        self.volume.nuclear_solid_target = 0

        # set fluid change rate
        self.fluid_change_rate = self.unlysed_fluid_change_rate


class NecrosisSwell(Phase):
    """
    Swelling part of the necrosis process
    """

    def __init__(self, index: int = 0, previous_phase_index: int = 0, next_phase_index: int = 1, dt: float = 0.1,
                 time_unit: str = "min", name: str = "Necrotic (swelling)", division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = None,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate: float = None,
                 nuclear_biomass_change_rate: float = None, calcification_rate: float = None,
                 relative_rupture_volume: float = None,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None):

        # default parameters

        # this phase will use by default a custom transition check, so this is here to avoid issues
        _phase_duration = 9e99

        # default volume parameters
        _cytoplasm_biomass_change_rate = 0.0032 / 60.0
        _nuclear_biomass_change_rate = 0.013 / 60.0
        _fluid_change_rate = 0.67 / 60.0
        _calcification_rate = 0.0042 / 60.0
        _relative_rupture_volume = 2

        if phase_duration is None:
            phase_duration = _phase_duration

        if cytoplasm_biomass_change_rate is None:
            cytoplasm_biomass_change_rate = _cytoplasm_biomass_change_rate

        if nuclear_biomass_change_rate is None:
            nuclear_biomass_change_rate = _nuclear_biomass_change_rate

        if fluid_change_rate is None:
            fluid_change_rate = _fluid_change_rate

        if calcification_rate is None:
            calcification_rate = _calcification_rate

        if relative_rupture_volume is None:
            relative_rupture_volume = _relative_rupture_volume

        if entry_function is None:
            entry_function = self._standard_necrosis_entry_function
            entry_function_args = [None]

        if transition_to_next_phase is None:
            transition_to_next_phase = self._necrosis_transition_function
            transition_to_next_phase_args = [None]

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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)

    def _standard_necrosis_entry_function(self, *none):

        # the cell wants to degrade the solids and swell by osmosis
        self.volume.target_fluid_fraction = 1
        self.volume.nuclear_solid_target = 0
        self.volume.cytoplasm_solid_target = 0

        self.volume.target_cytoplasm_to_nuclear_ratio = 0

        # set rupture volume

        self.volume.rupture_volume = self.volume.relative_rupture_volume * self.volume.total

    def _necrosis_transition_function(self, *none):
        return self.volume.total > self.volume.rupture_volume


class NecrosisLysed(Phase):
    """
    Ruptured necrotic cell
    """

    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 99, dt: float = 0.1,
                 time_unit: str = "min", name: str = "Necrotic (lysed)", division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = True, fixed_duration: bool = True, phase_duration: float = None,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 transition_to_next_phase=None, transition_to_next_phase_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_biomass_change_rate: float = None,
                 nuclear_biomass_change_rate: float = None, calcification_rate: float = None,
                 relative_rupture_volume: float = None,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None):

        # default parameters

        _phase_duration = 60 * 60 * 24  # 60 days, the cell should disappear naturally before then,
        # but if it hasn't we do it

        _cytoplasm_biomass_change_rate = 0.0032 / 60.0
        _nuclear_biomass_change_rate = 0.013 / 60.0
        _fluid_change_rate = 0.050 / 60.0
        _calcification_rate = 0.0042 / 60.0
        _relative_rupture_volume = 9e99

        if phase_duration is None:
            phase_duration = _phase_duration

        if cytoplasm_biomass_change_rate is None:
            cytoplasm_biomass_change_rate = _cytoplasm_biomass_change_rate

        if nuclear_biomass_change_rate is None:
            nuclear_biomass_change_rate = _nuclear_biomass_change_rate

        if fluid_change_rate is None:
            fluid_change_rate = _fluid_change_rate

        if calcification_rate is None:
            calcification_rate = _calcification_rate

        if relative_rupture_volume is None:
            relative_rupture_volume = _relative_rupture_volume

        if entry_function is None:
            entry_function = self._standard_lysis_entry_function
            entry_function_args = [None]

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
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume)

    def _standard_lysis_entry_function(self, *none):
        self.volume.target_fluid_fraction = 0
        self.volume.nuclear_solid_target = 0
        self.volume.cytoplasm_solid_target = 0

        self.volume.target_cytoplasm_to_nuclear_ratio = 0

        # set rupture volume

        self.volume.rupture_volume = self.volume.relative_rupture_volume * self.volume.total


if __name__ == '__main__':
    test_ki = Ki67Positive(dt=0.1)
    print(test_ki.index)
