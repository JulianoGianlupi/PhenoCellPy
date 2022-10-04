from .. import cycles


# from .. import phenotypes


def add_cycle_to_object(o: object, cycle: str or cycles.cycle, name: str = "unnamed",
                        dt: float = 1, time_unit: str = "min", phases: list = None, quiescent_phase= None):

    # todo: check if `o` can have custom attributes
    # if not hasattr(o, __dict__):
    #     raise

    if not isinstance(cycle, cycles.Cycle) and cycle not in cycles.cycle_names:
        raise ValueError("Expected `cycle` parameter to be either an initialized instance of Cycle class or the name "
                         "of a "
                         f"pre-defined cycle. Got {cycle}")

    if type(cycle) == str:
        cycle = cycles.get_cycle_by_name(cycle)
        cycle = cycle(name=name, dt=dt, time_unit=time_unit, phases=phases, quiescent_phase=quiescent_phase)

    setattr(o, "cycle", cycle)

