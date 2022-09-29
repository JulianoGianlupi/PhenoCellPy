import Phenotypes.phases as Phases

class Cycle:
    def __init__(self, name: str = "unnamed", dt: float = 1, time_unit: str = "min", phases: list = None,
                 quiescent_phase: Phases.Phase or bool = None):

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
            self.quiecent_phase = Phases.QuiescentPhase(dt=self.dt)
        elif quiescent_phase is not None and not quiescent_phase:
            self.quiecent_phase = False
        else:
            self.quiecent_phase = quiescent_phase
        self.current_phase = self.phases[0]
        self.time_in_cycle = 0

    def time_step_cycle(self):

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
        divides = self.current_phase.division_at_phase_exit
        dies = self.current_phase.removal_at_phase_exit
        self.set_phase(self.current_phase.next_phase_index)
        return True, dies, divides

    def set_phase(self, idx):
        self.current_phase = self.phases[idx]
        self.current_phase.time_in_phase = 0

    def go_to_quiescence(self):
        if self.quiecent_phase is not None and not self.quiecent_phase:
            return
        self.current_phase = self.quiecent_phase
        self.current_phase.time_in_phase = 0


class SimpleLiveCycle(Cycle):
    def __init__(self, time_unit: str = "min", name: str = "simple_live", dt=1):
        phases = [Phases.Phase(index=0, previous_phase_index=0, next_phase_index=0, time_unit=time_unit, name="alive",
                        division_at_phase_exit=True, phase_duration=60, dt=dt)]
        super().__init__(name=name, time_unit=time_unit, phases=phases, quiescent_phase=False, dt=dt)


class Ki67Basic(Cycle):
    def __init__(self, name="Ki67 Basic", dt=0.1):  # todo: init signature
        Ki67_positive = Phases.Phase()  # todo: define as a stand alone class
        Ki67_negative = Phases.Phase()  # todo: define as a stand alone class

        phases = [Ki67_negative, Ki67_positive]

        super().__init__(name=name, dt=dt, phases=phases, quiescent_phase=False)