"""
BSD 3-Clause License

Copyright (c) 2023, Juliano Ferrari Gianlupi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from numpy import exp
from numpy.random import uniform

from PhenoCellPy.cell_volume import CellVolumes


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
        senescence). See `time_step_phase`'s documentation for further explanation.

    transition_to_next_phase(*args)
        One of the default transition functions (`_transition_to_next_phase_deterministic`,
        `_transition_to_next_phase_stochastic`) or a user defined function. If user defined, it will be called with
        `check_transition_to_next_phase_function_args` as args. Must return a bool denoting if the transition occurs or not.

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
        Optional function that returns true if the cell should exit the cell cycle and enter senescence

    user_phase_time_step(*args)
        User-defined function to be executed with the time-step


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

    :param next_phase_index: Index of the phase proceeding this phase in the list of phases that forms a phenotype
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

    :param fixed_duration: Boolean setting the transition from this phase to the next to be deterministic (True) or
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
    and enter senescence
    :type arrest_function: function

    :param arrest_function_args: Args for `arrest_function`
    :type arrest_function_args: list or tuple

    :param check_transition_to_next_phase_function: Default or custom function that returns if the cell should
    advance to the next phase in the phenotype. If left as None the phase will pick either
    `_transition_to_next_phase_deterministic` or `_transition_to_next_phase_stochastic` depending on the value of
    `fixed_duration`. :type check_transition_to_next_phase_function: function

    :param check_transition_to_next_phase_function_args: Args for `transition_to_next_phase`
    :type check_transition_to_next_phase_function_args: list or tuple

    :param simulated_cell_volume: Volume of the simulated cell (e.g., a CompuCell3D or Tissue Forge cell)
    :type simulated_cell_volume: float

    :param cytoplasm_volume_change_rate: Change rate for the cytoplasmic volume. volume/`time_unit` units.
    `cytoplasm_volume_change_rate` >= 0. Passed to the `CellVolume` attribute class
    :type cytoplasm_volume_change_rate: float

    :param nuclear_volume_change_rate: Change rate for the nuclear volume. volume/`time_unit` units.
    `nuclear_volume_change_rate` >= 0. Passed to the `CellVolume` attribute class
    :type nuclear_volume_change_rate: float

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

    :param user_phase_time_step_args: args for `user_phase_time_step`
    :type user_phase_time_step_args: list or tuple

    Attributes
    ----------

    time_in_phase : float
        Time spent in this phase

    volume : class:cell_volume.CellVolumes
        Cell volume submodel

    """

    def __init__(self, index: int = None, previous_phase_index: int = None, next_phase_index: int = None,
                 dt: float = None, time_unit: str = "min", space_unit="micrometer", name: str = None,
                 division_at_phase_exit: bool = False, removal_at_phase_exit: bool = False,
                 fixed_duration: bool = False, phase_duration: float = 10, entry_function=None,
                 entry_function_args: list = None, exit_function=None, exit_function_args: list = None,
                 arrest_function=None, arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None, calcification_rate=None,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None,
                 relative_rupture_volume=None, user_phase_time_step=None, user_phase_time_step_args=(None,)):
        """

        :param space_unit:
        :type space_unit:
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
        and enter senescence
        :type arrest_function: function

        :param arrest_function_args: Args for `arrest_function`
        :type arrest_function_args: list or tuple

        :param check_transition_to_next_phase_function: Default or custom function that returns if the cell should
        advance to the next
        phase in the phenotype. If left as None the phase will pick either `_transition_to_next_phase_deterministic`
        or `_transition_to_next_phase_stochastic` depending on the value of `fixed_duration`.
        :type check_transition_to_next_phase_function: function

        :param check_transition_to_next_phase_function_args: Args for `transition_to_next_phase`
        :type check_transition_to_next_phase_function_args: list or tuple

        :param simulated_cell_volume: Volume of the simulated cell (e.g., a CompuCell3D or Tissue Forge cell)
        :type simulated_cell_volume: float

        :param cytoplasm_volume_change_rate: Change rate for the cytoplasmic volume. volume/`time_unit` units.
        `cytoplasm_volume_change_rate` >= 0. Passed to the `CellVolume` attribute class
        :type cytoplasm_volume_change_rate: float

        :param nuclear_volume_change_rate: Change rate for the nuclear volume. volume/`time_unit` units.
        `nuclear_volume_change_rate` >= 0. Passed to the `CellVolume` attribute class
        :type nuclear_volume_change_rate: float

        :param calcification_rate: Rate of calcification of the cell. volume/`time_unit` units. `calcification_rate`
        >= 0
        Passed to the `CellVolume` attribute class
        :type calcification_rate: float

        :param target_fluid_fraction: Fraction of the cell volume it will attempt to keep as fluid.
        0 <= `target_fluid_fraction` <= 1. Passed to the `CellVolume` attribute class
        :type target_fluid_fraction: float

        :param nuclear_fluid: Initial volume of the fluid part of the nucleus. Passed to the `CellVolume` attribute
        class
        :type nuclear_fluid: float

        :param nuclear_solid: Initial volume of the solid part of the nucleus. Passed to the `CellVolume` attribute
        class
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

        :param user_phase_time_step_args: args for `user_phase_time_step`
        :type user_phase_time_step_args: list or tuple

        """

        if index is None:
            self.index = 0  # int
        else:
            self.index = index

        self.previous_phase_index = previous_phase_index

        self.next_phase_index = next_phase_index

        self.time_unit = time_unit
        self.space_unit = space_unit

        if dt is None or dt <= 0:
            raise ValueError(f"'dt' must be greater than 0. Got {dt}.")
        self.dt = dt

        if name is None:
            self.name = "unnamed"  # string: phase's name
        else:
            self.name = name

        self.division_at_phase_exit = division_at_phase_exit  # bool flagging for division
        self.removal_at_phase_exit = removal_at_phase_exit  # bool flagging for removal (e.g., death)

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

        self.arrest_function = arrest_function  # function determining if cell will exit cell cycle and become senescent
        self.arrest_function_args = arrest_function_args

        if self.arrest_function is not None and type(self.arrest_function_args) != list:
            raise TypeError("Arrest function defined but no args given. Was expecting "
                            f"'arrest_function_args' to be a list, got {type(arrest_function_args)}.")

        if check_transition_to_next_phase_function is None:
            self.check_transition_to_next_phase_function_args = [None]
            if fixed_duration:
                self.check_transition_to_next_phase_function = self._check_transition_to_next_phase_deterministic
            else:
                self.check_transition_to_next_phase_function = self._check_transition_to_next_phase_stochastic
        else:
            if type(check_transition_to_next_phase_function_args) != list:
                raise TypeError("Custom exit function selected but no args given. Was expecting "
                                "'check_transition_to_next_phase_function_args' to be a list, got "
                                f"{type(check_transition_to_next_phase_function_args)}.")
            self.check_transition_to_next_phase_function_args = check_transition_to_next_phase_function_args
            self.check_transition_to_next_phase_function = check_transition_to_next_phase_function

        if simulated_cell_volume is None:
            self.simulated_cell_volume = 1
        else:
            self.simulated_cell_volume = simulated_cell_volume

        # the default rates are reference values for MCF-7, in 1/min
        if cytoplasm_volume_change_rate is None:
            self.cytoplasm_volume_change_rate = 0.27 / 60.0
        else:
            self.cytoplasm_volume_change_rate = cytoplasm_volume_change_rate

        if nuclear_volume_change_rate is None:
            self.nuclear_volume_change_rate = 0.33 / 60.0
        else:
            self.nuclear_volume_change_rate = nuclear_volume_change_rate
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

        self.user_phase_time_step = user_phase_time_step

        if self.user_phase_time_step is not None and \
                not (type(user_phase_time_step_args) == list or type(user_phase_time_step_args) == tuple):
            raise ValueError(
                f"`user_phase_time_step` is defined but `user_phase_time_step_args` is not list or "
                f"tuple.\nGot {type(user_phase_time_step_args)} instead")

        self.user_phase_time_step_args = user_phase_time_step_args

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
        self.volume.update_volume(self.dt, self.fluid_change_rate, self.nuclear_volume_change_rate,
                                  self.cytoplasm_volume_change_rate, self.calcification_rate)

    def _check_transition_to_next_phase_stochastic(self, *none):
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

    def _check_transition_to_next_phase_deterministic(self, *none):
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
        cycle (i.e., leaves the cycle; goes to senescence). If the cell doesn't senesce, this function checks if the
        cell
        should transition to the next phase.

        :return: tuple. First element of tuple: bool denoting if the cell moves to the next phase. Second element:
        denotes if the cell leaves the cell cycle and enters senescence.
        """
        self.time_in_phase += self.dt

        self.update_volume()

        if self.user_phase_time_step is not None:
            self.user_phase_time_step(*self.user_phase_time_step_args)

        if self.arrest_function is not None:
            exit_phenotype = self.arrest_function(*self.exit_function_args)
            go_to_next_phase_in_phenotype = False
            return go_to_next_phase_in_phenotype, exit_phenotype
        else:
            exit_phenotype = False

        go_to_next_phase_in_phenotype = self.check_transition_to_next_phase_function(*self.check_transition_to_next_phase_function_args)

        if go_to_next_phase_in_phenotype and self.exit_function is not None:
            self.exit_function(*self.exit_function_args)
            return go_to_next_phase_in_phenotype, exit_phenotype
        return go_to_next_phase_in_phenotype, exit_phenotype

    def _double_target_volume(self, *none):
        """

        Doubles the cell volume submodel (:class:`PhenoCellPy.cell_volume`) target volumes. Used by several cell cycle
        models to double the cell volume before mitosis

        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return: No return
        """
        self.volume.nuclear_solid_target *= 2
        self.volume.cytoplasm_solid_target *= 2

    def _halve_target_volume(self, *none):
        """

        Halves the cell volume submodel (:class:`PhenoCellPy.cell_volume`) target volumes. Used by several cell cycle
        models to halve the cell volume after mitosis

        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return: No return
        """
        self.volume.cytoplasm_solid_target /= 2
        self.volume.nuclear_solid_target /= 2

    def __str__(self):
        return f"{self.name} phase"


class SenescentPhase(Phase):
    """
    Default Senescent Phase. Inherits :class:`Phase`

    This senescent phase class is meant to be "outside" whatever phenotype progression is being used.

    """

    def __init__(self, index: int = 9999, previous_phase_index: int = None, next_phase_index: int = 9999,
                 dt: float = None, time_unit: str = "min", space_unit="micrometer", name: str = "senescent",
                 division_at_phase_exit: bool = False, removal_at_phase_exit: bool = False,
                 fixed_duration: bool = True, phase_duration: float = 60 * 24 * 60, entry_function=None,
                 entry_function_args: list = None, exit_function=None, exit_function_args: list = None,
                 arrest_function=None, arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate=0, nuclear_volume_change_rate=0, calcification_rate=0,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None,
                 relative_rupture_volume=None, user_phase_time_step=None, user_phase_time_step_args=None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit, name=name,
                         division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)
        return


class Ki67Negative(Phase):
    """
    Inherits :class:`Phase`. Defines Ki 67- quiescent phase.

    This is a quiescent phenotype for cells that are replicating. Ki67 is a protein marker associated with
    proliferation.
    Transition to the next phase is set to be stochastic (the phase does not use a fixed duration) by default. Default
    expected phase duration is 4.59h, the phase transition rate is, therefore, dt/4.59 1/h.
    This phase does not calcify the cell. The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939
    """

    def __init__(self, index: int = 0, previous_phase_index: int = 1, next_phase_index: int = 1, dt: float = 0.1,
                 time_unit: str = "min", space_unit="micrometer", name: str = "Ki 67-",
                 division_at_phase_exit: bool = False, removal_at_phase_exit: bool = False,
                 fixed_duration: bool = False, phase_duration: float = 4.59 * 60, entry_function=None,
                 entry_function_args: list = None, exit_function=None, exit_function_args: list = None,
                 arrest_function=None, arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None, calcification_rate=None,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None,
                 relative_rupture_volume=None, user_phase_time_step=None, user_phase_time_step_args=None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit,
                         name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)


class Ki67Positive(Phase):
    """

    Inherits :class:`Phase`. Defines Ki 67+ proliferating phase.

    This is a proliferating phenotype for cells that are replicating. Ki67 is a protein marker associated with
    proliferation. Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by
    default.
    Default phase duration is 15.5h. By default, if no user defined custom entry function is defined (i.e.,
    `entry_function=None`), this phase will set its entry function to be :class:`Phase._double_target_volume`.
    By default, will set the volume change rates to be [change in volume]/[phase duration]. This phase does not calcify
    the cell.

    The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939

    """

    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 0, dt: float = 0.1,
                 time_unit: str = "min", space_unit="micrometer", name: str = "Ki 67+",
                 division_at_phase_exit: bool = True, removal_at_phase_exit: bool = False, fixed_duration: bool = True,
                 phase_duration: float = 15.5 * 60.0, entry_function=None, entry_function_args: list = None,
                 exit_function=None, exit_function_args: list = None, arrest_function=None,
                 arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None, calcification_rate=None,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None,
                 relative_rupture_volume=None, user_phase_time_step=None, user_phase_time_step_args=None):

        if entry_function is None:
            entry_function = self._double_target_volume
            entry_function_args = [None]
        elif entry_function == False:  # CANNOT BE SIMPLIFIED TO `not entry_function`. Because not None == True
            pass
        elif type(entry_function_args) != list and type(entry_function_args) != tuple:
            raise TypeError("'entry_function' was defined but no valid value for 'entry_function_args' was given. "
                            "Expected "
                            f"list or tuple got {type(entry_function_args)}")

        if exit_function == False:  # CANNOT be changed to not exit_function!!! not None => True, None == False => False
            exit_function = None
            exit_function_args = [None]
        elif exit_function is None:
            exit_function = self._halve_target_volume
            exit_function_args = [None]
        elif type(exit_function_args) != list and type(exit_function_args) != tuple:
            raise TypeError("'exit_function' was defined but no  valid value for 'entry_function_args' was given. "
                            "Expected "
                            f"list or tuple got {type(exit_function_args)}")
        if cytoplasm_volume_change_rate is None and cytoplasm_fluid is not None and cytoplasm_solid is not None:
            cytoplasm_volume_change_rate = (cytoplasm_fluid + cytoplasm_solid) / (phase_duration / dt)

        elif cytoplasm_volume_change_rate is None and cytoplasm_fluid is not None:
            cytoplasm_volume_change_rate = cytoplasm_fluid / (phase_duration / dt)

        elif cytoplasm_volume_change_rate is None and cytoplasm_solid is not None:
            cytoplasm_volume_change_rate = cytoplasm_solid / (phase_duration / dt)

        elif cytoplasm_volume_change_rate is None:
            cytoplasm_volume_change_rate = 1
        else:
            cytoplasm_volume_change_rate = cytoplasm_volume_change_rate

        if nuclear_volume_change_rate is None and cytoplasm_fluid is not None and cytoplasm_solid is not None:
            nuclear_volume_change_rate = (nuclear_fluid + nuclear_solid) / (phase_duration / dt)

        elif nuclear_volume_change_rate is None and cytoplasm_fluid is not None:
            nuclear_volume_change_rate = nuclear_fluid / (phase_duration / dt)

        elif nuclear_volume_change_rate is None and cytoplasm_solid is not None:
            nuclear_volume_change_rate = nuclear_solid / (phase_duration / dt)

        elif nuclear_volume_change_rate is None:
            nuclear_volume_change_rate = 1
        else:
            nuclear_volume_change_rate = nuclear_volume_change_rate

        if fluid_change_rate is None and cytoplasm_fluid is not None and nuclear_fluid is not None:
            fluid_change_rate = (cytoplasm_fluid + nuclear_fluid) / (phase_duration / dt)
        elif fluid_change_rate is None and cytoplasm_fluid is not None:
            fluid_change_rate = cytoplasm_fluid / (phase_duration / dt)
        elif fluid_change_rate is None and nuclear_fluid is not None:
            fluid_change_rate = nuclear_fluid / (phase_duration / dt)
        else:
            fluid_change_rate = 1

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit,
                         name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)


class Ki67PositivePreMitotic(Ki67Positive):
    """

    Inherits :class:`Ki67Positive`. Defines Ki 67+ pre-mitotic proliferating phase. Only difference to
    :class:`Ki67Positive` is the phase length.

    This is a proliferating phenotype for cells that are replicating. Ki67 is a protein marker associated with
    proliferation. Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by
    default.
    Default phase duration is 13h. By default, if no user defined custom entry function is defined (i.e.,
    `entry_function=None`), this phase will set its entry function to be :class:`Phase._double_target_volume`

    The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939

    """

    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 2, dt: float = 0.1,
                 time_unit: str = "min", space_unit="micrometer", name: str = "Ki 67+ pre-mitotic",
                 division_at_phase_exit: bool = True, removal_at_phase_exit: bool = False, fixed_duration: bool = True,
                 phase_duration: float = 13.0 * 60.0, entry_function=None, entry_function_args: list = None,
                 exit_function=None, exit_function_args: list = None, arrest_function=None,
                 arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None, calcification_rate=None,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None,
                 relative_rupture_volume=None, user_phase_time_step=None, user_phase_time_step_args=None):
        if entry_function is None:
            # otherwise it will be defaulted to the halving target volume function by Ki67Positive
            entry_function = False

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit,
                         name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)


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
                 time_unit: str = "min", space_unit="micrometer", name: str = "Ki 67+ post-mitotic",
                 division_at_phase_exit: bool = True, removal_at_phase_exit: bool = False, fixed_duration: bool = True,
                 phase_duration: float = 2.5 * 60.0, entry_function=None, entry_function_args: list = None,
                 exit_function=None, exit_function_args: list = None, arrest_function=None,
                 arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None, calcification_rate=None,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None,
                 relative_rupture_volume=None, user_phase_time_step=None, user_phase_time_step_args=None):

        if entry_function is None:
            entry_function = self._standard_Ki67_positive_postmit_entry_function
            entry_function_args = [None]
        elif type(entry_function_args) != list:
            raise TypeError("'entry_function' was defined but no value for 'entry_function_args' was given. Expected "
                            f"list got {type(entry_function_args)}")

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit,
                         name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)

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
                 time_unit: str = "min", space_unit="micrometer", name: str = "G0/G1",
                 division_at_phase_exit: bool = False, removal_at_phase_exit: bool = False,
                 fixed_duration: bool = False, phase_duration: float = 5.15 * 60.0, entry_function=None,
                 entry_function_args: list = None, exit_function=None, exit_function_args: list = None,
                 arrest_function=None, arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None, calcification_rate=None,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None,
                 relative_rupture_volume=None, user_phase_time_step=None, user_phase_time_step_args=None):
        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit,
                         name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)


class S(Phase):
    """
    Inherits :class:`Phase`. Defines S phase, it more representative of the growth phase than the inter-growth rest.

    This is a growth phenotype for cells that are replicating. Transition to the next phase is set to be stochastic
    (the phase does not use a fixed duration) by default. Default expected phase duration is 8h, the phase transition
    rate is, therefore, dt/8 1/h. By default, will set the volume change rates to be
    [change in volume]/[phase duration]. This phase does not calcify the cell. Reference phase duration from
    https://www.ncbi.nlm.nih.gov/books/NBK9876/
    """

    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = 2, dt: float = 0.1,
                 time_unit: str = "min", space_unit="micrometer", name: str = "S", division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 8 * 60.0,
                 entry_function=None, entry_function_args: list = None, exit_function=None,
                 exit_function_args: list = None, arrest_function=None, arrest_function_args: list = None,
                 check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None,
                 simulated_cell_volume: float = None, cytoplasm_volume_change_rate=None,
                 nuclear_volume_change_rate=None, calcification_rate=None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None, relative_rupture_volume=None,
                 user_phase_time_step=None, user_phase_time_step_args=None):

        if entry_function is None:
            entry_function = self._double_target_volume
            entry_function_args = [None]
        elif type(entry_function_args) != list:
            raise TypeError("'entry_function' was defined but no value for 'entry_function_args' was given. Expected "
                            f"list got {type(entry_function_args)}")

        if cytoplasm_volume_change_rate is None and cytoplasm_fluid is not None and cytoplasm_solid is not None:
            cytoplasm_volume_change_rate = (cytoplasm_fluid + cytoplasm_solid) / (phase_duration / dt)

        elif cytoplasm_volume_change_rate is None and cytoplasm_fluid is not None:
            cytoplasm_volume_change_rate = cytoplasm_fluid / (phase_duration / dt)

        elif cytoplasm_volume_change_rate is None and cytoplasm_solid is not None:
            cytoplasm_volume_change_rate = cytoplasm_solid / (phase_duration / dt)

        elif cytoplasm_volume_change_rate is None:
            cytoplasm_volume_change_rate = 1
        else:
            cytoplasm_volume_change_rate = cytoplasm_volume_change_rate

        if nuclear_volume_change_rate is None and cytoplasm_fluid is not None and cytoplasm_solid is not None:
            nuclear_volume_change_rate = (nuclear_fluid + nuclear_solid) / (phase_duration / dt)

        elif nuclear_volume_change_rate is None and cytoplasm_fluid is not None:
            nuclear_volume_change_rate = nuclear_fluid / (phase_duration / dt)

        elif nuclear_volume_change_rate is None and cytoplasm_solid is not None:
            nuclear_volume_change_rate = nuclear_solid / (phase_duration / dt)

        elif nuclear_volume_change_rate is None:
            nuclear_volume_change_rate = 1
        else:
            nuclear_volume_change_rate = nuclear_volume_change_rate

        if fluid_change_rate is None and cytoplasm_fluid is not None and nuclear_fluid is not None:
            fluid_change_rate = (cytoplasm_fluid + nuclear_fluid) / (phase_duration / dt)
        elif fluid_change_rate is None and cytoplasm_fluid is not None:
            fluid_change_rate = cytoplasm_fluid / (phase_duration / dt)
        elif fluid_change_rate is None and nuclear_fluid is not None:
            fluid_change_rate = nuclear_fluid / (phase_duration / dt)
        else:
            fluid_change_rate = 1

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit,
                         name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)


class G2M(Phase):
    """
    Inherits :class:`Phase`. Defines G2M phase, it more representative of the mitosis phase than the growth.

    This is a growth phenotype for cells that are replicating. Transition to the next phase is set to be stochastic
    (the phase does not use a fixed duration) by default. Default expected phase duration is 5h, the phase transition
    rate is, therefore, dt/5 1/h. By default, will set the volume change rates to be
    [change in volume]/[phase duration]. This phase does not calcify the cell. Reference phase duration from
    https://www.ncbi.nlm.nih.gov/books/NBK9876/
    """

    def __init__(self, index: int = 2, previous_phase_index: int = 1, next_phase_index: int = 0, dt: float = 0.1,
                 time_unit: str = "min", space_unit="micrometer", name: str = "G2/M",
                 division_at_phase_exit: bool = True, removal_at_phase_exit: bool = False, fixed_duration: bool = False,
                 phase_duration: float = 5 * 60.0, entry_function=None, entry_function_args: list = None,
                 exit_function=None, exit_function_args: list = None, arrest_function=None,
                 arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None, calcification_rate=None,
                 target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, fluid_change_rate=None,
                 relative_rupture_volume=None, user_phase_time_step=None, user_phase_time_step_args=None):
        if entry_function is None:
            entry_function = self._double_target_volume
            entry_function_args = [None]
        elif type(entry_function_args) != list and type(entry_function_args) != tuple:
            raise TypeError("'entry_function' was defined but no valid value for 'entry_function_args' was given. "
                            "Expected "
                            f"list or tuple got {type(entry_function_args)}")
        if exit_function is None:
            exit_function = self._halve_target_volume
            exit_function_args = [None]
        elif type(exit_function_args) != list and type(exit_function_args) != tuple:
            raise TypeError("'exit_function' was defined but no  valid value for 'entry_function_args' was given. "
                            "Expected "
                            f"list or tuple got {type(exit_function_args)}")

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit,
                         name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)


class Apoptosis(Phase):
    """
    Inherits :class:`Phase`. Defines apoptotic phenotype phase.

    This is a dead cell phenotype phase, the cell will shrink itself and should be removed from the simulation when this
    phase ends. Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by
    default. Default phase duration is 8.6h. By default, if no custom user defined entry function is used (i.e.,
    `entry_function=None`), entry function is set to :class:`Apoptosis._standard_apoptosis_entry`.
    :class:`Apoptosis._standard_apoptosis_entry` sets all the cell target volumes from :class:`PhenoCellPy.cell_volume`
    to 0. The default mass change rates are `cytoplasm_volume_change_rate = 1/60` [volume/min],
    `nuclear_volume_change_rate = 0.35 / 60` [volume/min], `fluid_change_rate = 3 / 60`. This phase does not calcify
    the cell.
    """

    def __init__(self, index: int = 0, previous_phase_index: int = 0, next_phase_index: int = 0, dt: float = 0.1,
                 time_unit: str = "min", space_unit="micrometer", name: str = "Apoptosis",
                 division_at_phase_exit: bool = False, removal_at_phase_exit: bool = True, fixed_duration: bool = True,
                 phase_duration: float = 8.6 * 60.0, entry_function=None, entry_function_args: list = None,
                 exit_function=None, exit_function_args: list = None, arrest_function=None,
                 arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate: float = 1 / 60, nuclear_volume_change_rate: float = 0.35 / 60,
                 calcification_rate: float = 0, relative_rupture_volume: float = 2, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=3 / 60, user_phase_time_step=None,
                 user_phase_time_step_args=None):

        if entry_function is None:
            entry_function = self._standard_apoptosis_entry
            entry_function_args = [None]

        if fluid_change_rate is not None:
            self.fluid_change_rate = fluid_change_rate
        else:
            self.fluid_change_rate = 3 / 60

        if relative_rupture_volume is None:
            self.relative_rupture_volume = 2
        else:
            self.relative_rupture_volume = relative_rupture_volume

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit,
                         name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)

    def _standard_apoptosis_entry(self, *none):
        """
        Zeroes all the cell's target volumes. Keeps the nuclear to cytoplasm ratio the same.

        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return:
        """

        # shrink cell
        self.volume.target_fluid_fraction = 0
        self.volume.cytoplasm_solid_target = 0
        self.volume.nuclear_solid_target = 0


class NecrosisSwell(Phase):
    """
    Inherits :class:`Phase`. Swelling part of the necrosis process.

    Represents the osmotic swell a necrotic cell goes through. By default, this phase uses a custom transition function
    (i.e., `check_transition_to_next_phase_functions=None`), it can be overwritten by a user defined one. The custom
    transition
    function is :class:`NecrosisSwell._necrosis_transition_function`, it returns true when the cell becomes bigger than
    its rupture volume. The default relative rupture volume is 2, i.e., the cell ruptures after doubling in volume.
    By default, if no custom user defined entry function is used (i.e., `entry_function=None`), entry function is set
    to :class:`NecrosisSwell._standard_necrosis_entry_function`. It zeroes the solid target volumes and the target
    cytoplasm to nuclear ratio, and sets the target fluid fraction to 1. This causes the cell to increase its volume.
    The default volume change rates are `cytoplasm_volume_change_rate = 0.0032 / 60.0`,
    `nuclear_volume_change_rate = 0.013 / 60.0`, `fluid_change_rate = 0.67 / 60.0`,
    `calcification_rate = 0.0042 / 60.0`. This phase does calcify the cell.
    """

    def __init__(self, index: int = 0, previous_phase_index: int = 0, next_phase_index: int = 1, dt: float = 0.1,
                 time_unit: str = "min", space_unit="micrometer", name: str = "Necrotic (swelling)",
                 division_at_phase_exit: bool = False, removal_at_phase_exit: bool = False,
                 fixed_duration: bool = False, phase_duration: float = None, entry_function=None,
                 entry_function_args: list = None, exit_function=None, exit_function_args: list = None,
                 arrest_function=None, arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate: float = None, nuclear_volume_change_rate: float = None,
                 calcification_rate: float = None, relative_rupture_volume: float = None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None, user_phase_time_step=None,
                 user_phase_time_step_args=None):

        # default parameters

        # this phase will use by default a custom transition check, so this is here to avoid issues
        _phase_duration = 9e99

        # default volume parameters
        _cytoplasm_volume_change_rate = 0.0032 / 60.0
        _nuclear_volume_change_rate = 0.013 / 60.0
        _fluid_change_rate = 0.67 / 60.0
        _calcification_rate = 0.0042 / 60.0
        _relative_rupture_volume = 2

        if phase_duration is None:
            phase_duration = _phase_duration

        if cytoplasm_volume_change_rate is None:
            cytoplasm_volume_change_rate = _cytoplasm_volume_change_rate

        if nuclear_volume_change_rate is None:
            nuclear_volume_change_rate = _nuclear_volume_change_rate

        if fluid_change_rate is None:
            fluid_change_rate = _fluid_change_rate

        if calcification_rate is None:
            calcification_rate = _calcification_rate

        if relative_rupture_volume is None:
            relative_rupture_volume = _relative_rupture_volume

        if entry_function is None:
            entry_function = self._standard_necrosis_entry_function
            entry_function_args = [None]

        if check_transition_to_next_phase_function is None:
            check_transition_to_next_phase_function = self._necrosis_transition_function
            check_transition_to_next_phase_function_args = [None]

        super().__init__(index=index, previous_phase_index=previous_phase_index, next_phase_index=next_phase_index,
                         dt=dt, time_unit=time_unit, space_unit=space_unit,
                         name=name, division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)

    def _standard_necrosis_entry_function(self, *none):
        """
        Responsible for causing the osmotic swell.

        Zeroes the solid target volumes, sets the target cytoplasm to nuclear
        ratio to 0, sets the target fluid fraction to 1, sets the rupture volume to be double the current total volume.
        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return: No return
        """

        # the cell wants to degrade the solids and swell by osmosis
        self.volume.target_fluid_fraction = 1
        self.volume.nuclear_solid_target = 0
        self.volume.cytoplasm_solid_target = 0

        self.volume.target_cytoplasm_to_nuclear_ratio = 0

        # set rupture volume

        self.volume.rupture_volume = self.volume.relative_rupture_volume * self.volume.total

    def _necrosis_transition_function(self, *none):
        """
        Custom phase transition function. The simulated cell should only change phase once it bursts (i.e., when its
        volume is above the rupture volume), it cares not how long or how little time it takes to reach that state.
        :param none: Not used. This is a custom transition function, therefore it has to have args
        :return: Flag for phase transition
        :rtype: bool
        """
        return self.volume.total > self.volume.rupture_volume


class NecrosisLysed(Phase):
    """
    Inherits :class:`Phase`. Ruptured necrotic cell

    Represents the already ruptured necrotic cell. The modeler should find a way to represent this fragmentary state,
    either create several tiny cells from the cell that ruptured, or create a disconnected cell.  By default, the
    transition to the next phase is deterministic and the phase lasts 60 days. The simulated cell should shrink and
    disappear before then, this is a safeguard to remove the cell "by hand" if it hasn't. By default, if no custom user
    defined entry function is used (i.e., `entry_function=None`), entry function is set to
    :class:`NecrosisLysed._standard_lysis_entry_function`. It zeroes all target volumes from
    :class:`PhenoCellPy.cell_volume`. The default volume change rates are:
    `cytoplasm_volume_change_rate = 0.0032 / 60.0`, `nuclear_volume_change_rate = 0.013 / 60.0`,
    `fluid_change_rate = 0.050 / 60.0`, `calcification_rate = 0.0042 / 60.0`. This phase calcifies the cell.

    """

    def __init__(self, index: int = 1, previous_phase_index: int = 0, next_phase_index: int = -1, dt: float = 0.1,
                 time_unit: str = "min", space_unit="micrometer", name: str = "Necrotic (lysed)",
                 division_at_phase_exit: bool = False, removal_at_phase_exit: bool = True, fixed_duration: bool = True,
                 phase_duration: float = None, entry_function=None, entry_function_args: list = None,
                 exit_function=None, exit_function_args: list = None, arrest_function=None,
                 arrest_function_args: list = None, check_transition_to_next_phase_function=None,
                 check_transition_to_next_phase_function_args: list = None, simulated_cell_volume: float = None,
                 cytoplasm_volume_change_rate: float = None, nuclear_volume_change_rate: float = None,
                 calcification_rate: float = None, relative_rupture_volume: float = None, target_fluid_fraction=None,
                 nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                 cytoplasm_solid=None, cytoplasm_solid_target=None, target_cytoplasm_to_nuclear_ratio=None,
                 calcified_fraction=None, fluid_change_rate=None, user_phase_time_step=None,
                 user_phase_time_step_args=None):

        # default parameters

        _phase_duration = 60 * 60 * 24  # 60 days, the cell should disappear naturally before then,
        # but if it hasn't we do it

        _cytoplasm_volume_change_rate = 0.0032 / 60.0
        _nuclear_volume_change_rate = 0.013 / 60.0
        _fluid_change_rate = 0.050 / 60.0
        _calcification_rate = 0.0042 / 60.0
        _relative_rupture_volume = 9e99

        if phase_duration is None:
            phase_duration = _phase_duration

        if cytoplasm_volume_change_rate is None:
            cytoplasm_volume_change_rate = _cytoplasm_volume_change_rate

        if nuclear_volume_change_rate is None:
            nuclear_volume_change_rate = _nuclear_volume_change_rate

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
                         dt=dt, time_unit=time_unit, space_unit=space_unit, name=name,
                         division_at_phase_exit=division_at_phase_exit,
                         removal_at_phase_exit=removal_at_phase_exit, fixed_duration=fixed_duration,
                         phase_duration=phase_duration, entry_function=entry_function,
                         entry_function_args=entry_function_args, exit_function=exit_function,
                         exit_function_args=exit_function_args, arrest_function=arrest_function,
                         arrest_function_args=arrest_function_args,
                         check_transition_to_next_phase_function=check_transition_to_next_phase_function,
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                         nuclear_volume_change_rate=nuclear_volume_change_rate, calcification_rate=calcification_rate,
                         target_fluid_fraction=target_fluid_fraction, nuclear_fluid=nuclear_fluid,
                         nuclear_solid=nuclear_solid, nuclear_solid_target=nuclear_solid_target,
                         cytoplasm_fluid=cytoplasm_fluid, cytoplasm_solid=cytoplasm_solid,
                         cytoplasm_solid_target=cytoplasm_solid_target,
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
                         calcified_fraction=calcified_fraction, fluid_change_rate=fluid_change_rate,
                         relative_rupture_volume=relative_rupture_volume, user_phase_time_step=user_phase_time_step,
                         user_phase_time_step_args=user_phase_time_step_args)

    def _standard_lysis_entry_function(self, *none):
        """
        Zeroes all the cell's target volumes. Also zeroes the nuclear to cytoplasm ratio.
        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return:
        """
        self.volume.target_fluid_fraction = 0
        self.volume.nuclear_solid_target = 0
        self.volume.cytoplasm_solid_target = 0

        self.volume.target_cytoplasm_to_nuclear_ratio = 0

        # set rupture volume

        self.volume.rupture_volume = self.volume.relative_rupture_volume * self.volume.total


if __name__ == '__main__':
    dt = 1
    phase = Phase(dt=dt)
    sen = SenescentPhase(dt=dt)
    ki67n = Ki67Negative(dt=dt)
    test_ki67p = Ki67Positive(dt=dt)
    ki67ppre = Ki67PositivePreMitotic(dt=dt)
    ki67ppos = Ki67PositivePostMitotic(dt=dt)
    g0g1 = G0G1(dt=dt)
    s = S(dt=dt)
    g2m = G2M(dt=dt)
    ap = Apoptosis(dt=dt)
    necsw = NecrosisSwell(dt=dt)
    necLys = NecrosisLysed(dt=dt)


    def grow_phase_transition(*args):
        return args[0] >= args[1] and args[2] > args[4]


    def double_target_volumes(self, *none):
        self.volume.nuclear_solid_target *= 2
        self.volume.cytoplasm_solid_target *= 2

    custom = Phase(index=1, previous_phase_index=0, next_phase_index=2, dt=dt,
                   time_unit="min", space_unit="micrometer", name="custom",
                   division_at_phase_exit=False, removal_at_phase_exit=False, fixed_duration=True,
                   phase_duration=120, entry_function=double_target_volumes,
                   entry_function_args=[None],
                   exit_function=None, arrest_function=None,
                   check_transition_to_next_phase_function=grow_phase_transition,
                   check_transition_to_next_phase_function_args=[0, 9, 0, 9],
                   simulated_cell_volume=1,
                   cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None,
                   calcification_rate=None, target_fluid_fraction=None, nuclear_fluid=None,
                   nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                   cytoplasm_solid=None, cytoplasm_solid_target=None,
                   target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None,
                   fluid_change_rate=None, relative_rupture_volume=None,
                   user_phase_time_step=None, user_phase_time_step_args=(None,))

    # print(test_ki67p.index)
    for _ in range(1000):
        # print(phase.name, phase.time_step_phase(), phase.volume.total)
        # print(qui.name, qui.time_step_phase(), qui.volume.total)
        # print(ki67n.name, ki67n.time_step_phase(), ki67n.volume.total)
        # print(test_ki67p.name, test_ki67p.time_step_phase(), test_ki67p.volume.total)
        # print(ki67ppre.name, ki67ppre.time_step_phase(), ki67ppre.volume.total)
        # print(ki67ppos.name, ki67ppos.time_step_phase(), ki67ppos.volume.total)
        # print(g0g1.name, g0g1.time_step_phase(), g0g1.volume.total)
        # print(s.name, s.time_step_phase(), s.volume.total)
        # print(g2m.name, g2m.time_step_phase(), g2m.volume.total)
        # print(ap.name, ap.time_step_phase(), ap.volume.total)
        # print(necsw.name, necsw.time_step_phase(), necsw.volume.total)
        # print(necLys.name, necLys.time_step_phase(), necLys.volume.total)
        print(custom.name, custom.time_step_phase(), custom.volume.total)
