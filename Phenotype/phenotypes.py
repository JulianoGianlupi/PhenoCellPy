
class Phase:
    def __init__(self, index: int = None, previous_phase_index: int = None, next_phase_index: int = None,
                 name: str = None, division_at_phase_exit: bool = False, removal_at_phase_exit: bool = False,
                 entry_function=None, exit_function=None, arrest_function=None):
        if index is None:
            self.index = 0  # int
        else:
            self.index = index

        self.previous_phase_index = previous_phase_index

        self.next_phase_index = next_phase_index

        if name is None:
            self.name = "unnamed"  # string: phase's name
        else:
            self.name = name

        self.division_at_phase_exit = division_at_phase_exit  # bool flagging for division
        self.removal_at_phase_exit = removal_at_phase_exit  # bool flagging for removal (death?)

        self.entry_function = entry_function  # function to be executed upon entering this phase

        self.exit_function = exit_function  # function to be executed just before exiting this phase

        self.arrest_function = arrest_function  # function determining if
