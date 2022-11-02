

class BetterCellVolumes:  # todo: make it work for what I need
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

        self.__ff = (self.__cf + self.__nf)/(self.__nf+self.__ns+self.__cf+self.__cs)
        self.__nff = self.__nf / self.nuclear
        self.__cff = self.__cf / self.cytoplasm

    @property
    def total(self):
        return self.cytoplasm + self.nuclear

    @property
    def fluid(self):
        return self.cytoplasm_fluid + self.nuclear_fluid

    @fluid.setter
    def fluid(self, value):
        value = value if value >= 0 else 0
        self.nuclear_fluid = value * self.nuclear/self.total
        self.cytoplasm_fluid = value * self.cytoplasm/self.total

    @property
    def solid(self):
        return self.cytoplasm_solid + self.nuclear_solid

    @property
    def fluid_fraction(self):
        return self.__ff

    @fluid_fraction.setter
    def fluid_fraction(self, value):
        if 0 <= value <= 1:
            self.__ff = value
        elif value < 0:
            self.__ff = 0
        else:
            self.__ff = 1

    @property
    def nuclear_fluid_fraction(self):
        return self.__nff

    @nuclear_fluid_fraction.setter
    def nuclear_fluid_fraction(self, value):
        if 0 <= value <= 1:
            self.__nff = value
        elif value < 0:
            self.__nff = 0
        else:
            self.__nff = 1

    @property
    def cytoplasm_fluid_fraction(self):
        return self.__cff

    @cytoplasm_fluid_fraction.setter
    def cytoplasm_fluid_fraction(self, value):
        if 0 <= value <= 1:
            self.__cff = value
        elif value < 0:
            self.__cff = 0
        else:
            self.__cff = 1

    @property
    def solid_fraction(self):
        return self.solid / self.total

    @property
    def cytoplasm_to_nuclear_ratio(self):
        return self.cytoplasm / self.nuclear

    @property
    def cytoplasm(self):
        return self.cytoplasm_fluid + self.cytoplasm_solid

    @cytoplasm.setter
    def cytoplasm(self, value):
        value = value if value >= 0 else 0
        self.cytoplasm_fluid = self.cytoplasm_fluid_fraction * value
        self.cytoplasm_solid = (1 - self.cytoplasm_fluid_fraction) * value

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

    @nuclear.setter
    def nuclear(self, value):
        value = value if value >= 0 else 0
        self.nuclear_fluid = self.nuclear_fluid_fraction * value
        self.nuclear_solid = (1 - self.nuclear_fluid_fraction) * value

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

    # def update_fluid

    def update_volume(self, dt, cytoplasm_change_rate, nuclear_change_rate, calcification_rate):
        self.update_cytoplasm(dt, cytoplasm_change_rate)
        self.update_nuclear(dt, nuclear_change_rate)
        self.update_calcified(dt, calcification_rate)


class CellVolumes:
    def __init__(self, target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None):

        if target_fluid_fraction is None:
            self.target_fluid_fraction = 1
        else:
            self.target_fluid_fraction = target_fluid_fraction

        if nuclear_fluid is None:
            self.nuclear_fluid = 0
        else:
            self.nuclear_fluid = nuclear_fluid

        if nuclear_solid is None:
            self.nuclear_solid = 0
        else:
            self.nuclear_solid = nuclear_solid

        if nuclear_solid_target is None:
            self.nuclear_solid_target = 0
        else:
            self.nuclear_solid_target = nuclear_solid_target

        if cytoplasm_fluid is None:
            self.cytoplasm_fluid = 1
        else:
            self.cytoplasm_fluid = cytoplasm_fluid

        if cytoplasm_solid is None:
            self.cytoplasm_solid = 0
        else:
            self.cytoplasm_solid = cytoplasm_solid

        if cytoplasm_solid_target is None:
            self.cytoplasm_solid_target = 0
        else:
            self.cytoplasm_solid_target = cytoplasm_solid_target

        if target_cytoplasm_to_nuclear_ratio is None:
            self.target_cytoplasm_to_nuclear_ratio = 0
        else:
            self.target_cytoplasm_to_nuclear_ratio = target_cytoplasm_to_nuclear_ratio

        if calcified_fraction is None:
            self.calcified_fraction = 0
        else:
            self.calcified_fraction = calcified_fraction

        self.fluid = self.cytoplasm_fluid + self.nuclear_fluid

        self.cytoplasm = self.cytoplasm_fluid + self.cytoplasm_solid

        self.solid = self.cytoplasm_solid + self.nuclear_solid

        self.nuclear = self.nuclear_fluid + self.nuclear_solid

        self.total = self.nuclear + self.cytoplasm


    @property
    def fluid(self):
        return self.__fluid

    @fluid.setter
    def fluid(self, value):
        self.__fluid = value if value > 0 else 0

    @property
    def target_fluid_fraction(self):
        return self.__tff

    @target_fluid_fraction.setter
    def target_fluid_fraction(self, value):
        if value < 0:
            value = 0
        elif value > 1:
            value = 1
        self.__tff = value

    @property
    def total(self):
        return self.__total

    @total.setter
    def total(self, value):
        self.__total = value if value >= 0 else 0

    @property
    def nuclear_fluid(self):
        return self.__nuclear_fluid

    @nuclear_fluid.setter
    def nuclear_fluid(self, value):
        self.__nuclear_fluid = value if value >= 0 else 0

    @property
    def nuclear(self):
        return self.__nuclear

    @nuclear.setter
    def nuclear(self, value):
        self.__nuclear = value if value >= 0 else 0

    @property
    def cytoplasm_fluid(self):
        return self.__cytoplasm_fluid

    @cytoplasm_fluid.setter
    def cytoplasm_fluid(self, value):
        self.__cytoplasm_fluid = value if value >= 0 else 0

    @property
    def nuclear_solid(self):
        return self.__nuclear_solid

    @nuclear_solid.setter
    def nuclear_solid(self, value):
        self.__nuclear_solid = value if value >= 0 else 0

    @property
    def nuclear_solid_target(self):
        return self.__nst

    @nuclear_solid_target.setter
    def nuclear_solid_target(self, value):
        self.__nst = value if value >= 0 else 0

    @property
    def cytoplasm_solid_target(self):
        return self.__cst

    @cytoplasm_solid_target.setter
    def cytoplasm_solid_target(self, value):
        self.__cst = value if value >= 0 else 0

    @property
    def target_cytoplasm_to_nuclear_ratio(self):
        return self.__tctnr

    @target_cytoplasm_to_nuclear_ratio.setter
    def target_cytoplasm_to_nuclear_ratio(self, value):
        self.__tctnr = value if value >= 0 else 0

    @property
    def cytoplasm_solid(self):
        return self.__cytoplasm_solid

    @cytoplasm_solid.setter
    def cytoplasm_solid(self, value):
        self.__cytoplasm_solid = value if value >= 0 else 0

    @property
    def solid(self):
        return self.__solid

    @solid.setter
    def solid(self, value):
        self.__solid = value if value >= 0 else 0

    @property
    def cytoplasm(self):
        return self.__cytoplasm

    @cytoplasm.setter
    def cytoplasm(self, value):
        self.__cytoplasm = value if value >= 0 else 0

    @property
    def calcified_fraction(self):
        return self.__calc_frac

    @calcified_fraction.setter
    def calcified_fraction(self, value):
        if value > 1:
            value = 1
        elif value < 0:
            value = 0
        self.__calc_frac = value

    @property
    def fluid_fraction(self):
        return self.__fluid_fraction

    @fluid_fraction.setter
    def fluid_fraction(self, value):
        self.__fluid_fraction = value if value >= 0 else 0

    def update_volume(self, dt, fluid_change_rate, nuclear_biomass_change_rate, cytoplasm_biomass_change_rate,
                      calcification_rate):
        self.fluid += dt * fluid_change_rate * (self.target_fluid_fraction * self.total - self.fluid)

        self.nuclear_fluid = (self.nuclear / (self.total + 1e-12)) * self.fluid

        self.cytoplasm_fluid = self.fluid - self.nuclear_fluid

        self.nuclear_solid += dt * nuclear_biomass_change_rate * (self.nuclear_solid_target - self.nuclear_solid)

        self.cytoplasm_solid_target = self.target_cytoplasm_to_nuclear_ratio * self.nuclear_solid_target

        self.cytoplasm_solid += dt * cytoplasm_biomass_change_rate * (self.cytoplasm_solid_target -
                                                                      self.cytoplasm_solid)

        self.solid = self.nuclear_solid + self.cytoplasm_solid  # maybe this could be a pure property?

        self.nuclear = self.nuclear_solid + self.nuclear_fluid

        self.cytoplasm = self.cytoplasm_fluid + self.cytoplasm_solid

        self.calcified_fraction = dt * calcification_rate * (1 - self.calcified_fraction)

        self.total = self.cytoplasm + self.nuclear

        self.fluid_fraction = self.fluid / (self.total + 1e-12)


