
class Phase:
    """
    Generic cell cycle phase class.

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
    :param exit_function:
    :param arrest_function:
    """
    def __init__(self, index: int = None, previous_phase_index: int = None, next_phase_index: int = None,
                 time_unit: str = "min", name: str = None, division_at_phase_exit: bool = False,
                 removal_at_phase_exit: bool = False, fixed_duration: bool = False, phase_duration: float = 10,
                 entry_function=None, exit_function=None, arrest_function=None):

        if index is None:
            self.index = 0  # int
        else:
            self.index = index

        self.previous_phase_index = previous_phase_index

        self.next_phase_index = next_phase_index

        self.time_unit = time_unit

        if name is None:
            self.name = "unnamed"  # string: phase's name
        else:
            self.name = name

        self.division_at_phase_exit = division_at_phase_exit  # bool flagging for division
        self.removal_at_phase_exit = removal_at_phase_exit  # bool flagging for removal (death?)

        self.fixed_duration = fixed_duration

        self.phase_duration = phase_duration

        self.entry_function = entry_function  # function to be executed upon entering this phase

        self.exit_function = exit_function  # function to be executed just before exiting this phase

        self.arrest_function = arrest_function  # function determining if cell will exit cell cycle and become quiescent


if __name__ == '__main__':
    pass

    # testPhase = Phase(fixed_duration="True")
