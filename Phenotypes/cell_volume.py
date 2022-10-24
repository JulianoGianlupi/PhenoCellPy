class CellVolumes:
    def __init__(self, cytoplasm=None, target_cytoplasm=None, target_cytoplasm_fluid_fraction=None,
                 target_nuclear=None, target_nuclear_fluid_fraction=None, nuclear=None, calcified_fraction=None):

        if target_cytoplasm is None:
            self.__tc = .5
        else:
            if target_cytoplasm <= 0:
                raise ValueError(f"`target_cytoplasm` must be > 0, got {target_cytoplasm}")
            self.__tc = target_cytoplasm

        if target_cytoplasm_fluid_fraction is None:
            self.__tcff = 1
        else:
            if not 0 <= target_cytoplasm_fluid_fraction <= 1:
                raise ValueError(f"`target_cytoplasm_fluid_fraction` must be in range [0,1], got "
                                 f"{target_cytoplasm_fluid_fraction}")
            self.__tcff = target_cytoplasm_fluid_fraction

        if cytoplasm is None:
            self.__cf = self.__tcff * self.__tc
            self.__cs = (1 - self.__tcff) * self.__tc
        else:
            if cytoplasm < 0:
                raise ValueError(f"`cytoplasm` must be >= 0, got {cytoplasm}")
            self.__cf = self.__tcff * cytoplasm
            self.__cs = (1 - self.__tcff) * cytoplasm

        if target_nuclear is None:
            self.__tn = .5
        else:
            if target_nuclear < 0:
                raise ValueError(f"`target_nuclear` must be >= 0, got {target_nuclear}")
            self.__tn = target_nuclear

        if target_nuclear_fluid_fraction is None:
            self.__tnff = 1
        else:
            if not 0 <= target_nuclear_fluid_fraction <= 1:
                raise ValueError(f"`target_nuclear_fluid_fraction` must be in range [0,1], got "
                                 f"{target_nuclear_fluid_fraction}")
            self.__tnff = target_nuclear_fluid_fraction

        if nuclear is None:
            self.__nf = self.__tnff * self.__tc
            self.__ns = (1 - self.__tnff) * self.__tc
        else:
            if nuclear < 0:
                raise ValueError(f"`nuclear` must be >= 0, got {nuclear}")
            self.__nf = self.__tnff * nuclear
            self.__ns = (1 - self.__tnff) * nuclear

        if calcified_fraction is None:
            self.__cal_frac = 0
        else:
            if not 0 <= calcified_fraction <= 1:
                raise ValueError(f"`calcified_fraction` must be in range [0,1], got {calcified_fraction}")
            self.__cal_frac = calcified_fraction

    @property
    def total(self):
        return self.cytoplasm + self.nuclear

    @property
    def fluid(self):
        return self.cytoplasm_fluid + self.nuclear_fluid

    @property
    def solid(self):
        return self.cytoplasm_solid + self.nuclear_solid

    @property
    def fluid_fraction(self):
        return self.fluid / self.total

    @property
    def solid_fraction(self):
        return self.solid / self.total

    @property
    def cytoplasm_to_nuclear_ratio(self):
        return self.cytoplasm / self.nuclear

    @property
    def cytoplasm(self):
        return self.cytoplasm_fluid + self.cytoplasm_solid

    @property
    def target_cytoplasm(self):
        return self.__tc

    @target_cytoplasm.setter
    def target_cytoplasm(self, value):
        self.__tc = value if value >= 0 else 0

    @property
    def cytoplasm_fluid(self):
        return self.__cf

    @cytoplasm_fluid.setter
    def cytoplasm_fluid(self, value):
        self.__cf = value if value >= 0 else 0

    @property
    def cytoplasm_solid(self):
        return self.__cs

    @cytoplasm_solid.setter
    def cytoplasm_solid(self, value):
        self.__cs = value if value >= 0 else 0

    @property
    def target_cytoplasm_fluid_fraction(self):
        return self.__tcff

    @target_cytoplasm_fluid_fraction.setter
    def target_cytoplasm_fluid_fraction(self, value):
        if 0 <= value <= 1:
            self.__tcff = value
        elif value < 0:
            self.__tcff = 0
        else:
            self.__tcff = 1

    @property
    def nuclear(self):
        return self.nuclear_fluid + self.nuclear_solid

    @property
    def target_nuclear(self):
        return self.__tn

    @target_nuclear.setter
    def target_nuclear(self, value):
        self.__tn = value if value >= 0 else 0

    @property
    def nuclear_fluid(self):
        return self.__nf

    @nuclear_fluid.setter
    def nuclear_fluid(self, value):
        self.__nf = value if value >= 0 else 0

    @property
    def nuclear_solid(self):
        return self.__ns

    @nuclear_solid.setter
    def nuclear_solid(self, value):
        self.__ns = value if value >= 0 else 0

    @property
    def target_nuclear_fluid_fraction(self):
        return self.__tnff

    @target_nuclear_fluid_fraction.setter
    def target_nuclear_fluid_fraction(self, value):
        if 0 <= value <= 1:
            self.__tnff = value
        elif value < 0:
            self.__tnff = 0
        else:
            self.__tnff = 1

    @property
    def calcified_fraction(self):
        return self.__cal_frac

    @calcified_fraction.setter
    def calcified_fraction(self, value):
        if 0 <= value <= 1:
            self.__cal_frac = value
        elif value < 0:
            self.__cal_frac = 0
        else:
            self.__cal_frac = 1

    def update_cytoplasm(self, dt, change_rate):
        self.cytoplasm_fluid += dt * change_rate * (self.target_cytoplasm_fluid_fraction * self.target_cytoplasm -
                                                    self.cytoplasm_fluid)
        self.cytoplasm_solid += dt * change_rate * ((1 - self.target_cytoplasm_fluid_fraction) * self.target_cytoplasm -
                                                    self.cytoplasm_solid)

    def update_nuclear(self, dt, change_rate):
        self.nuclear_fluid += dt * change_rate * (self.target_nuclear_fluid_fraction * self.target_nuclear -
                                                  self.nuclear_fluid)
        self.nuclear_solid += dt * change_rate * ((1 - self.target_nuclear_fluid_fraction) * self.target_nuclear -
                                                  self.nuclear_solid)

    def update_calcified(self, dt, change_rate):
        self.calcified_fraction = dt * change_rate * (1 - self.calcified_fraction)

    def update_volume(self, dt, cytoplasm_change_rate, nuclear_change_rate, calcification_rate):
        self.update_cytoplasm(dt, cytoplasm_change_rate)
        self.update_nuclear(dt, nuclear_change_rate)
        self.update_calcified(dt, calcification_rate)
