

class Phase:
    def __init__(self, index=None, name=None, division_at_phase_exit=False, removal_at_phase_exit=False,
                 entry_function=None):
        if index is None:
            self.index = 0
        else:
            self.index = index
        if name is None:
            self.name = "unnamed"
        else:
            self.name = name

        self.division_at_phase_exit = division_at_phase_exit
        self.removal_at_phase_exit = removal_at_phase_exit

        self.entry_function = entry_function


