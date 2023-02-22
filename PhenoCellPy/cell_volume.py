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
from scipy.integrate import odeint
from numpy import array
from copy import deepcopy

class CellVolumes:
    """
    Cell volume class, evolves the cell volume and its subvolumes

    Methods:
    --------

    update_volume(dt, fluid_change_rate, nuclear_volume_change_rate, cytoplasm_volume_change_rate, calcification_rate)
        Time steps the volume model by dt. Updates all the cell's subvolumes based on the set targets by using scipy's
        odeint and the general ODE form of the volume relaxation `volume_relaxation`

    Static Methods:
    ---------------
    volume_relaxation(current_volume, t, rate, target_volume)
        General form of the volume dynamics ODE. Passed to scipy's odeint to be solved.


    Parameters:
    -----------

    :param target_fluid_fraction: Fraction of the cell volume that should be liquid
    :type target_fluid_fraction: float in [0, 1]
    :param nuclear_solid_target: How much of the nuclear volume should be solid
    :type nuclear_solid_target: float
    :param cytoplasm_solid_target: How much of the cytoplasm volume should be solid
    :type cytoplasm_solid_target: float
    :param target_cytoplasm_to_nuclear_ratio: How big the ratio "cytoplasm volume / nuclear volume" should be
    :type target_cytoplasm_to_nuclear_ratio: float
    :param relative_rupture_volume: Relative volume at which the cell should burst, set at the start. To do anything
    should be checked by Phase or Phenotype
    :type relative_rupture_volume: float

    Attributes:
    -----------

    :param nuclear_fluid: How much of the nuclear volume is fluid
    :type nuclear_fluid: float
    :param nuclear_solid: How much of the nuclear volume is solid
    :type nuclear_solid: float
    :param cytoplasm_fluid: How much of the cytoplasm volume is fluid
    :type cytoplasm_fluid: float
    :param cytoplasm_solid: How much of the cytoplasm volume is solid
    :type cytoplasm_solid: float
    :param calcified_fraction: How much of the cell is calcified
    :type calcified_fraction: float in range [0, 1]

    Default parameters:
    ------------------

    The defaults values below are reference parameter values for MCF-7, in cubic microns
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-7

    _total = 2494
    _fluid_fraction = .75
    _fluid = _fluid_fraction * _total
    _solid = _total - _fluid
    _nuclear = 540
    _nuclear_fluid = _fluid_fraction * _nuclear
    _nuclear_solid = _nuclear - _nuclear_fluid
    _cytoplasm = _total - _nuclear
    _cytoplasm_fluid = _fluid_fraction * _cytoplasm
    _cytoplasm_solid = _cytoplasm - _cytoplasm_fluid
    _calcified_fraction = 0
    _relative_rupture_volume = 100

    """

    def __init__(self, target_fluid_fraction=None, nuclear_fluid=None, nuclear_solid=None, nuclear_solid_target=None,
                 cytoplasm_fluid=None, cytoplasm_solid=None, cytoplasm_solid_target=None,
                 target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None, relative_rupture_volume=None,
                 time_unit="min", space_unit="micrometer"):
        """

        :param time_unit:
        :type time_unit:
        :param space_unit:
        :type space_unit:
        :param target_fluid_fraction: Fraction of the cell volume that should be liquid
        :type target_fluid_fraction: float in [0, 1]
        :param nuclear_fluid: How much of the nuclear volume is fluid
        :type nuclear_fluid: float
        :param nuclear_solid: How much of the nuclear volume is solid
        :type nuclear_solid: float
        :param nuclear_solid_target: How much of the nuclear volume should be solid
        :type nuclear_solid_target: float
        :param cytoplasm_fluid: How much of the cytoplasm volume is fluid
        :type cytoplasm_fluid: float
        :param cytoplasm_solid: How much of the cytoplasm volume is solid
        :type cytoplasm_solid: float
        :param cytoplasm_solid_target: How much of the cytoplasm volume should be solid
        :type cytoplasm_solid_target: float
        :param target_cytoplasm_to_nuclear_ratio: How big the ratio "cytoplasm volume / nuclear volume" should be
        :type target_cytoplasm_to_nuclear_ratio: float
        :param calcified_fraction: How much of the cell is calcified
        :type calcified_fraction: float in range [0, 1]
        :param relative_rupture_volume: Relative volume at which the cell should burst, set at the start. To do anything
        should be checked by Phase or Phenotype
        :type relative_rupture_volume: float
        """
        # The defaults values below are reference parameter values for MCF-7, in cubic
        # https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-7
        _total = 2494
        _fluid_fraction = .75
        _fluid = _fluid_fraction * _total
        _solid = _total - _fluid
        _nuclear = 540
        _nuclear_fluid = _fluid_fraction * _nuclear
        _nuclear_solid = _nuclear - _nuclear_fluid
        _cytoplasm = _total - _nuclear
        _cytoplasm_fluid = _fluid_fraction * _cytoplasm
        _cytoplasm_solid = _cytoplasm - _cytoplasm_fluid
        _calcified_fraction = 0
        _relative_rupture_volume = 100

        # setting class parameters

        self.time_unit = time_unit
        self.space_unit = space_unit

        if target_fluid_fraction is None:
            self.target_fluid_fraction = _fluid_fraction
        else:
            if not 0 <= target_fluid_fraction <= 1:
                raise ValueError(f"`target_fluid_fraction` must be in range [0, 1]. Got {target_fluid_fraction}")
            self.target_fluid_fraction = target_fluid_fraction

        if nuclear_fluid is None:
            self.nuclear_fluid = _nuclear * self.target_fluid_fraction
        else:
            if nuclear_fluid < 0:
                raise ValueError(f"`nuclear_fluid` must be >=0. Got {nuclear_fluid}")
            self.nuclear_fluid = nuclear_fluid

        if nuclear_solid is None:
            self.nuclear_solid = _nuclear * (1 - self.target_fluid_fraction)
        else:
            if nuclear_solid < 0:
                raise ValueError(f"`nuclear_solid` must be >=0. Got {nuclear_solid}")
            self.nuclear_solid = nuclear_solid

        if nuclear_solid_target is None:
            self.nuclear_solid_target = self.nuclear_solid
        else:
            if nuclear_solid_target < 0:
                raise ValueError(f"`nuclear_solid_target` must be >=0. Got {nuclear_solid_target}")
            self.nuclear_solid_target = nuclear_solid_target

        if cytoplasm_fluid is None:
            self.cytoplasm_fluid = _cytoplasm * self.target_fluid_fraction
        else:
            self.cytoplasm_fluid = cytoplasm_fluid

        if cytoplasm_solid is None:
            self.cytoplasm_solid = _cytoplasm * (1 - self.target_fluid_fraction)
        else:
            self.cytoplasm_solid = cytoplasm_solid

        if cytoplasm_solid_target is None:
            self.cytoplasm_solid_target = self.cytoplasm_solid
        else:
            self.cytoplasm_solid_target = cytoplasm_solid_target

        self.cytoplasm = self.cytoplasm_fluid + self.cytoplasm_solid
        self.nuclear = self.nuclear_fluid + self.nuclear_solid

        if target_cytoplasm_to_nuclear_ratio is None:
            self.target_cytoplasm_to_nuclear_ratio = self.cytoplasm / (1e-16 + self.nuclear)
        else:
            self.target_cytoplasm_to_nuclear_ratio = target_cytoplasm_to_nuclear_ratio

        if calcified_fraction is None:
            self.calcified_fraction = _calcified_fraction
        else:
            self.calcified_fraction = calcified_fraction

        if relative_rupture_volume is None:
            self.relative_rupture_volume = _relative_rupture_volume
        else:
            self.relative_rupture_volume = relative_rupture_volume

        self.fluid = self.cytoplasm_fluid + self.nuclear_fluid

        self.solid = self.cytoplasm_solid + self.nuclear_solid

        self.total = self.nuclear + self.cytoplasm

        self.fluid_fraction = self.fluid / self.total

        self.rupture_volume = self.relative_rupture_volume * self.total

    @staticmethod
    def volume_relaxation(current_volume, t, rate, target_volume):
        """

        :param current_volume: Current volume value to be relaxed
        :type current_volume: float
        :param t: time points to solve the ODE at. Always [0, dt]
        :type t: numpy array
        :param rate: volume change rate for this volume
        :type rate: float
        :param target_volume: target volume to relax towards
        :type target_volume: float
        :return dvdt: ODE computed at times of t
        :rtype: float
        """
        dvdt = rate * (target_volume - current_volume)
        return dvdt

    def update_volume(self, dt, fluid_change_rate, nuclear_volume_change_rate, cytoplasm_volume_change_rate,
                      calcification_rate):
        """
        Updates the cell volume attributes.

        Dynamic volumes (class parameters) updated as first order differential equations with a set rate and set target.
        Other volumes (class attributes) are set as ratios/relations of the dynamic volumes.

        The dynamic volumes change is, unless stated otherwise,
        new_volume = volume + dt * change rate * (target -  volume)
        and are updated by a forward Euler method

        First the total fluid is updated with rate `fluid_change_rate` and target `target_fluid_fraction`.

        Then the nuclear and cytoplasm fluid volumes are set

        The next dynamic volume update is the solid nuclear with rate `fluid_change_rate` and target
        `target_fluid_fraction`.

        The solid target is then updated as `target_cytoplasm_to_nuclear_ratio` * `nuclear_solid_target`

        The cytoplasm solid is dynamically updated with rate `cytoplasm_volume_change_rate` and target
        `cytoplasm_solid_target`

        The solid, nuclear, and cytoplasm volumes are set

        The calcified fraction is dynamically calculated with rate `calcification_rate` and no rate

        Finally the total volume and fluid fractions are set

        :param dt: Time step length
        :param fluid_change_rate: How fast can the cell change the fluid fraction of its volume
        :param nuclear_volume_change_rate: How fast can the cell change its nuclear volume
        :param cytoplasm_volume_change_rate: How fast can the cell change its cytoplasm volume
        :param calcification_rate: The rate at which the cell is calcifying
        :return: None
        """

        dt_array = array([0, dt])

        # self.fluid += dt * fluid_change_rate * (self.target_fluid_fraction * self.total - self.fluid)
        # (current_volume, t, rate, target_volume)
        self.fluid = odeint(self.volume_relaxation, self.fluid, dt_array,
                            args=(fluid_change_rate, self.target_fluid_fraction * self.total))[-1][0]

        self.nuclear_fluid = (self.nuclear / (self.total + 1e-12)) * self.fluid

        self.cytoplasm_fluid = self.fluid - self.nuclear_fluid

        # self.nuclear_solid += dt * nuclear_volume_change_rate * (self.nuclear_solid_target - self.nuclear_solid)
        self.nuclear_solid = odeint(self.volume_relaxation, self.nuclear_solid, dt_array,
                                    args=(nuclear_volume_change_rate, self.nuclear_solid_target))[-1][0]

        self.cytoplasm_solid_target = self.target_cytoplasm_to_nuclear_ratio * self.nuclear_solid_target

        # self.cytoplasm_solid += dt * cytoplasm_volume_change_rate * (self.cytoplasm_solid_target -
        #                                                               self.cytoplasm_solid)
        # (current_volume, t, rate, target_volume)
        self.cytoplasm_solid = odeint(self.volume_relaxation, self.cytoplasm_solid, dt_array,
                                      args=(cytoplasm_volume_change_rate, self.cytoplasm_solid_target))[-1][0]

        self.solid = self.nuclear_solid + self.cytoplasm_solid  # maybe this could be a pure property?

        self.nuclear = self.nuclear_solid + self.nuclear_fluid

        self.cytoplasm = self.cytoplasm_fluid + self.cytoplasm_solid

        # self.calcified_fraction += dt * calcification_rate * (1 - self.calcified_fraction)
        # (current_volume, t, rate, target_volume)
        self.calcified_fraction = odeint(self.volume_relaxation, self.calcified_fraction, dt_array,
                                         args=(calcification_rate, 1))[-1][0]

        self.total = self.cytoplasm + self.nuclear

        self.fluid_fraction = self.fluid / (self.total + 1e-12)

    def copy(self):
        return deepcopy(self)

    @property
    def cytoplasm_to_nuclear_ratio(self):
        """Get the current cytoplasm to nuclear ratio"""
        return self.cytoplasm / (1e-16 + self.nuclear)

    @property
    def fluid(self):
        """Get the total fluid part of the cell volume"""
        return self._fluid

    @fluid.setter
    def fluid(self, value):
        """Negative volumes are not physical"""
        self._fluid = value if value > 0 else 0

    @property
    def target_fluid_fraction(self):
        """Get the current target fluid section"""
        return self._tff

    @target_fluid_fraction.setter
    def target_fluid_fraction(self, value):
        """A fraction of something must be in [0, 1]"""
        if value < 0:
            value = 0
        elif value > 1:
            value = 1
        self._tff = value

    @property
    def total(self):
        """Gets the total volume"""
        return self._total

    @total.setter
    def total(self, value):
        self._total = value if value >= 0 else 0

    @property
    def total_target(self):
        return (self.cytoplasm_solid_target + self.nuclear_solid_target)/(1 - self.target_fluid_fraction)

    @property
    def nuclear_fluid(self):
        """Gets the nuclear fluid volume"""
        return self._nuclear_fluid

    @nuclear_fluid.setter
    def nuclear_fluid(self, value):
        self._nuclear_fluid = value if value >= 0 else 0

    @property
    def nuclear(self):
        """Gets the nuclear volume"""
        return self._nuclear

    @nuclear.setter
    def nuclear(self, value):
        self._nuclear = value if value >= 0 else 0

    @property
    def cytoplasm_fluid(self):
        """Gets the fluid volume of the cytoplasm"""
        return self._cytoplasm_fluid

    @cytoplasm_fluid.setter
    def cytoplasm_fluid(self, value):
        self._cytoplasm_fluid = value if value >= 0 else 0

    @property
    def nuclear_solid(self):
        """Gets the solid volume of the nucleus"""
        return self._nuclear_solid

    @nuclear_solid.setter
    def nuclear_solid(self, value):
        self._nuclear_solid = value if value >= 0 else 0

    @property
    def nuclear_solid_target(self):
        """Gets the target volume for the solid part of the nucleus"""
        return self._nst

    @nuclear_solid_target.setter
    def nuclear_solid_target(self, value):
        self._nst = value if value >= 0 else 0

    @property
    def cytoplasm_solid_target(self):
        """Gets the target volume for the solid part of the cytoplasm"""
        return self._cst

    @cytoplasm_solid_target.setter
    def cytoplasm_solid_target(self, value):
        self._cst = value if value >= 0 else 0

    @property
    def target_cytoplasm_to_nuclear_ratio(self):
        return self._tctnr

    @target_cytoplasm_to_nuclear_ratio.setter
    def target_cytoplasm_to_nuclear_ratio(self, value):
        """Gets the target ratio cytoplasm/nucleus"""
        self._tctnr = value if value >= 0 else 0

    @property
    def cytoplasm_solid(self):
        """Gets the solid part of the cytoplasm volume"""
        return self._cytoplasm_solid

    @cytoplasm_solid.setter
    def cytoplasm_solid(self, value):
        self._cytoplasm_solid = value if value >= 0 else 0

    @property
    def solid(self):
        """Gets the solid volume of the cell"""
        return self._solid

    @solid.setter
    def solid(self, value):
        self._solid = value if value >= 0 else 0

    @property
    def cytoplasm(self):
        """Gets the cytoplasm volume"""
        return self._cytoplasm

    @cytoplasm.setter
    def cytoplasm(self, value):
        self._cytoplasm = value if value >= 0 else 0

    @property
    def calcified_fraction(self):
        """Gets the calcified fraction of the cell"""
        return self._calc_frac

    @calcified_fraction.setter
    def calcified_fraction(self, value):
        if value > 1:
            value = 1
        elif value < 0:
            value = 0
        self._calc_frac = value

    @property
    def fluid_fraction(self):
        """Gets how much of the cell volume if fluid as a fraction of the total"""
        return self._fluid_fraction

    @fluid_fraction.setter
    def fluid_fraction(self, value):
        self._fluid_fraction = value if value >= 0 else 0



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


if __name__ == "__main__":
    from numpy.random import uniform
    cv = CellVolumes()

    for i in range(300):
        cv.update_volume(.500, 2, 2, 2, 2)

        print(cv.total, (cv.nuclear_solid_target + cv.cytoplasm_solid_target) / (1 - cv.target_fluid_fraction))
        if i == 100:
            # cv.target_cytoplasm_to_nuclear_ratio *= 2
            # cv.target_fluid_fraction *= uniform(0.8, 1.2)
            cv.nuclear_solid_target *= 2
            cv.cytoplasm_solid_target *= 2
        if i == 200:
            # cv.target_cytoplasm_to_nuclear_ratio *= 2
            # cv.target_fluid_fraction *= uniform(0.8, 1.2)
            cv.nuclear_solid_target /= 2
            cv.cytoplasm_solid_target /= 2
