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
        if self.current_phase.entry_function is not None:
            self.current_phase.entry_function(self.current_phase.entry_function_args)
        self.current_phase.time_in_phase = 0
        return True, dies, divides

    def set_phase(self, idx):
        volume = self.current_phase.volume
        tg_volume = self.current_phase.target_volume
        self.current_phase = self.phases[idx]
        self.current_phase.volume = volume
        self.current_phase.target_volume = tg_volume
        self.current_phase.time_in_phase = 0

    def go_to_quiescence(self):
        if self.quiecent_phase is not None and not self.quiecent_phase:
            return
        volume = self.current_phase.volume
        self.quiecent_phase.volume = volume
        self.quiecent_phase.target_volume = volume
        self.current_phase = self.quiecent_phase
        self.current_phase.time_in_phase = 0


class SimpleLiveCycle(Cycle):
    def __init__(self, time_unit: str = "min", name: str = "Simple Live", dt=1):
        phases = [Phases.Phase(index=0, previous_phase_index=0, next_phase_index=0, dt=dt, time_unit=time_unit,
                               name="alive", division_at_phase_exit=True, phase_duration=60)]
        super().__init__(name=name, time_unit=time_unit, phases=phases, quiescent_phase=False, dt=dt)


class Ki67Basic(Cycle):
    def __init__(self, name="Ki67 Basic", dt=0.1, quiescent_phase=False, time_unit="min",
                 target_volumes: list = None, volumes: list = None, update_volume_rate=None):

        if target_volumes is None:
            target_volumes = [1, 1]
            volumes = [1, 1]

        Ki67_positive = Phases.Ki67Negative(index=1, dt=dt, previous_phase_index=0, next_phase_index=0,
                                            target_volume=target_volumes[1], volume=volumes[1], time_unit=time_unit,
                                            update_volume_rate=update_volume_rate)
        Ki67_negative = Phases.Ki67Positive(index=0, dt=dt, previous_phase_index=1, next_phase_index=1,
                                            target_volume=target_volumes[0], volume=volumes[0], time_unit=time_unit)

        phases = [Ki67_negative, Ki67_positive]

        super().__init__(name=name, dt=dt, phases=phases, quiescent_phase=quiescent_phase, time_unit=time_unit)


cycle_names = ["Simple Live", "Ki67 Basic"]


def get_cycle_by_name(name):

    if name not in cycle_names:
        raise ValueError(f"{name} is not a pre-defined cycle")

    if name == "Simple Live":
        return SimpleLiveCycle
    elif name == "Ki67 Basic":
        return Ki67Basic

    return Cycle


if __name__=="__main__":
    print(cycle_names)



