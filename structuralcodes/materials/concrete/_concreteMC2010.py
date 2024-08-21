"""The concrete class for Model Code 2020 Concrete Material."""

import typing as t
import warnings

from structuralcodes.codes import mc2010

from ._concrete import Concrete


class ConcreteMC2010(Concrete):
    """Concrete implementation for MC 2010."""

    # computed values
    _fcm: t.Optional[float] = None
    _fctm: t.Optional[float] = None
    _Eci: t.Optional[float] = None
    _fctkmin: t.Optional[float] = None
    _fctkmax: t.Optional[float] = None
    _Gf: t.Optional[float] = None
    _Eci: t.Optional[float] = None
    _alpha_cc: t.Optional[float] = None
    _eps_c1: t.Optional[float] = None
    _eps_cu1: t.Optional[float] = None
    _k_sargin: t.Optional[float] = None
    _eps_c2: t.Optional[float] = None
    _eps_cu2: t.Optional[float] = None
    _n_parabolic_rectangular: t.Optional[float] = None
    _eps_c3: t.Optional[float] = None
    _eps_cu3: t.Optional[float] = None

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400.0,
        gamma_c: t.Optional[float] = None,
        existing: bool = False,
        alpha_cc: t.Optional[float] = None,
        **kwargs,
    ):
        """Initializes a new instance of Concrete for MC 2010.

        Arguments:
            fck (float): Characteristic strength in MPa if concrete is not
                existing.

        Keyword Arguments:
            name (Optional(str)): A descriptive name for concrete.
            density (float): Density of material in kg/m3 (default: 2400).
            gamma_c (Optional(float)): The partial factor for concrete.
            existing (bool): The material is of an existing structure
                (default: False).
            alpha_cc (float, optional): A factor for considering long-term
                effects on the strength, and effects that arise from the way
                the load is applied.
        """
        del kwargs
        if name is None:
            name = f'C{round(fck):d}'
        super().__init__(
            fck=fck,
            name=name,
            density=density,
            existing=existing,
            gamma_c=gamma_c,
        )
        self._alpha_cc = alpha_cc

    def _reset_attributes(self):
        self._fcm = None
        self._fctm = None
        self._fctkmin = None
        self._fctkmax = None
        self._Gf = None
        self._Eci = None
        self._eps_c1 = None
        self._eps_cu1 = None
        self._k_sargin = None
        self._eps_c2 = None
        self._eps_cu2 = None
        self._n_parabolic_rectangular = None
        self._eps_c3 = None
        self._eps_cu3 = None

    @property
    def fcm(self) -> float:
        """Returns fcm in MPa.

        Returns:
            float: The mean compressive strength in MPa.
        """
        self._fcm = self._fcm or mc2010.fcm(self._fck)
        return self._fcm

    @fcm.setter
    def fcm(self, value: float):
        """Sets a user defined value for fcm.

        Arguments:
            value (float): The value of fcm in MPa.

        Raises:
            ValueError: If value is lower than fck.
        """
        if abs(value) <= self._fck:
            raise ValueError(
                (
                    'Mean compressive strength cannot be lower than',
                    'characteristic strength.\n',
                    'Current characteristing strength: ',
                    f'fck = {self._fck}.',
                    f'Current value: value = {value}',
                )
            )
        self._fcm = abs(value)

    @property
    def Eci(self) -> float:
        """Returns the modulus of elasticity in MPa at the concrete age of 28
        days.

        It is assumed a normal concrete with quartzite aggregates (alfa_e = 1)
        """
        self._Eci = self._Eci or mc2010.Eci(self.fcm)
        return self._Eci

    @Eci.setter
    def Eci(self, value: float):
        """Sets a user defined value for modulus of elasticity at the concrete
        age of 28 days, Eci.

        Arguments:
            value (float): The value of Eci in MPa.
        """
        if value < 1e4 or value > 1e5:
            warnings.warn(
                'A suspect value of Eci has been input.\n'
                'Please check Eci that should be in MPa ({value} given).'
            )
        self._Eci = abs(value)

    @property
    def fctm(self) -> float:
        """Returns fctm in MPa.

        Returns:
            float: The mean tensile strength in MPa.
        """
        self._fctm = self._fctm or mc2010.fctm(self._fck)
        return self._fctm

    @fctm.setter
    def fctm(self, value: float):
        """Sets a user defined value for fctm.

        Arguments:
            value (float): The value of fctm in MPa.
        """
        if value > 0.5 * self._fck:
            warnings.warn(
                'A suspect value of fctm has been input. Please check.'
            )
        self._fctm = abs(value)

    @property
    def fctkmin(self) -> float:
        """Returns fctkmin in MPa.

        Returns:
            float: The lower bound tensile strength in MPa.
        """
        self._fctkmin = self._fctkmin or mc2010.fctkmin(self.fctm)
        return self._fctkmin

    @fctkmin.setter
    def fctkmin(self, value: float):
        """Sets a user defined value for fctkmin.

        Arguments:
            value (float): The value of fctkmin in MPa.
        """
        self._fctkmin = abs(value)

    @property
    def fctkmax(self) -> float:
        """Returns fctkmax in MPa.

        Returns:
            float: The upper bound tensile strength in MPa.
        """
        self._fctkmax = self._fctkmax or mc2010.fctkmax(self.fctm)
        return self._fctkmax

    @fctkmax.setter
    def fctkmax(self, value: float):
        """Sets a user defined value for fctkmax.

        Arguments:
            value (float): The value of fctkmax in MPa.
        """
        self._fctkmax = abs(value)

    @property
    def Gf(self) -> float:
        """Fracture energy of concrete.

        Returns:
            float: The fracture energy in N/m.
        """
        self._Gf = self._Gf or mc2010.Gf(self._fck)
        return self._Gf

    @Gf.setter
    def Gf(self, value: float):
        """Sets a user defined value for fracture energy Gf.

        Arguments:
            value (float): The value of Gf in N/m.
        """
        self._Gf = abs(value)

    @property
    def gamma_c(self) -> float:
        """The partial factor for concrete."""
        return self._gamma_c or 1.5

    def fcd(self) -> float:
        """Return the design compressive strength in MPa.

        Returns:
            float: The design compressive strength of concrete in MPa.
        """
        # This method should perhaps become a property, but is left as a method
        # for now, to be consistent with other concretes.
        return mc2010.fcd(
            self.fck, alpha_cc=self.alpha_cc, gamma_c=self.gamma_c
        )

    @property
    def alpha_cc(self) -> float:
        """The alpha_cc factor."""
        # Here we should implement the interaction with the globally set
        # national annex. For now, we simply return the default value.
        return self._alpha_cc or 1.0

    @property
    def eps_c1(self) -> float:
        """Returns the strain at maximum compressive strength of concrete (fcm)
        for the Sargin constitutive law.

        Returns:
            float: The strain at maximum compressive strength of concrete.
        """
        self._eps_c1 = self._eps_c1 or mc2010.eps_c1(self._fck)
        return self._eps_c1

    @eps_c1.setter
    def eps_c1(self, value: float):
        """Sets a user defined value for strain at peak strenght for Sargin
        constitutive law.

        Arguments:
            value (float): The new value for eps_c1, no units.
        """
        if abs(value) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_c1 that should be a pure'
                ' number without units. Plase check ({value} given).'
            )
        self._eps_c1 = value

    @property
    def eps_cu1(self) -> float:
        """Returns the strain at concrete failure of concrete.

        Returns:
            float: The maximum strength at failure of concrete.
        """
        self._eps_cu1 = self._eps_cu1 or mc2010.eps_cu1(self._fck)
        return self._eps_cu1

    @eps_cu1.setter
    def eps_cu1(self, value: float):
        """Sets the nominal ultimate strain for Sargin constitutive law.

        Arguments:
            value (float): The new value for eps_cu1, no units.
        """
        if abs(value) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_cu1 that should be a pure'
                ' number without units. Plase check ({value} given).'
            )
        self._eps_cu1 = value

    @property
    def k_sargin(self) -> float:
        """Returns the coefficient for Sargin constitutive law.

        Returns:
            float: The plastic coefficient for Sargin law.
        """
        self._k_sargin = self._k_sargin or mc2010.k_sargin(self._fck)
        return self._k_sargin

    @k_sargin.setter
    def k_sargin(self, value: float):
        """Sets the the coefficient for Sargin constitutive law.

        Arguments:
            value (float): The new value for k, no units.

        Raises:
            ValueError: If value < 0.
        """
        if value < 0:
            raise ValueError(f'n should be a positive value ({value} given)')
        self._k_sargin = value

    @property
    def eps_c2(self) -> float:
        """Returns the strain at maximum compressive strength of concrete (fcd)
        for the Parabola-rectangle constitutive law.

        Returns:
            float: The strain at maximum compressive strength of concrete.
        """
        self._eps_c2 = self._eps_c2 or mc2010.eps_c2(self.fck)
        return self._eps_c2

    @eps_c2.setter
    def eps_c2(self, value: float):
        """Sets the strain at maximum compressive strength of concrete (fcd)
        for the Parabola-rectangle constitutive law.

        Arguments:
            value (float): The new value for eps_c2, no units.
        """
        if abs(value) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_c2 that should be a pure'
                ' number without units. Plase check ({value} given).'
            )
        self._eps_c2 = value

    @property
    def eps_cu2(self) -> float:
        """Returns the strain at concrete failure of concrete for the
        Parabola-rectangle constitutive law.

        Returns:
            float: The maximum strain at failure of concrete.
        """
        self._eps_cu2 = self._eps_cu2 or mc2010.eps_cu2(self.fck)
        return self._eps_cu2

    @eps_cu2.setter
    def eps_cu2(self, value: float):
        """Sets the strain at concrete failure of concrete for the
        Parabola-rectangle constitutive law.

        Arguments:
            value (float): The new value for eps_cu2, no units.
        """
        if abs(value) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_cu2 that should be a pure'
                ' number without units. Plase check ({value} given).'
            )
        self._eps_cu2 = value

    @property
    def n_parabolic_rectangular(self) -> float:
        """Returns the coefficient for Parabola-rectangle constitutive law.

        Returns:
            float: The exponent for Parabola-recangle law.
        """
        self._n_parabolic_rectangular = (
            self._n_parabolic_rectangular
            or mc2010.n_parabolic_rectangular(self.fck)
        )
        return self._n_parabolic_rectangular

    @n_parabolic_rectangular.setter
    def n_parabolic_rectangular(self, value: float):
        """Sets the coefficient for Parabola-rectangle constitutive law.

        Arguments:
            value (float): The new value for n, no units.

        Raises:
            ValueError: If value < 0.
        """
        if value < 0:
            raise ValueError(f'n should be a positive value ({value} given)')
        if value >= 5:
            warnings.warn(
                'A suspect value is input for eps_cu2 that should be a pure'
                ' number without units. Plase check ({value} given).'
            )
        self._n_parabolic_rectangular = value

    @property
    def eps_c3(self) -> float:
        """Returns the strain at maximum compressive strength of concrete (fcd)
        for the Bi-linear constitutive law.

        Returns:
            float: The strain at maximum compressive strength of concrete.
        """
        self._eps_c3 = self._eps_c3 or mc2010.eps_c3(self.fck)
        return self._eps_c3

    @eps_c3.setter
    def eps_c3(self, value: float):
        """Sets the strain at maximum compressive strength of concrete (fcd)
        for the Bi-linear constitutive law.

        Arguments:
            value (float): The new value for eps_c3, no units.
        """
        if abs(value) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_c3 that should be a pure'
                ' number without units. Plase check ({value} given).'
            )
        self._eps_c3 = value

    @property
    def eps_cu3(self) -> float:
        """Returns the strain at concrete failure of concrete for the Bi-linear
        constitutive law.

        Returns:
            float: The maximum strain at failure of concrete.
        """
        self._eps_cu3 = self._eps_cu3 or mc2010.eps_cu3(self.fck)
        return self._eps_cu3

    @eps_cu3.setter
    def eps_cu3(self, value: float):
        """Sets the strain at concrete failure of concrete for the Bi-linear
        constitutive law.

        Arguments:
            value (float): The new value for eps_cu3, no units.
        """
        if abs(value) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_cu3 that should be a pure'
                ' number without units. Plase check ({value} given).'
            )
        self._eps_cu3 = value

    def __elastic__(self) -> dict:
        """Returns kwargs for creating an elastic constitutive law."""
        return {'E': self.Eci}

    def __bilinearcompression__(self) -> dict:
        """Returns kwargs for Bi-linear constitutive law."""
        return {
            'fc': self.fcd(),
            'eps_c': self.eps_c3,
            'eps_cu': self.eps_cu3,
        }

    def __parabolarectangle__(self) -> dict:
        """Returns kwargs for creating a parabola rectangle const law."""
        return {
            'fc': self.fcd(),
            'eps_0': self.eps_c2,
            'eps_u': self.eps_cu2,
            'n': self.n_parabolic_rectangular,
        }

    def __sargin__(self) -> dict:
        """Returns kwargs for creating a Sargin const law."""
        return {
            'fc': self.fcd(),
            'eps_c1': self.eps_c1,
            'eps_cu1': self.eps_cu1,
            'k': self.k_sargin,
        }

    def __popovics__(self) -> dict:
        """Returns kwargs for creating a Sargin const law."""
        return {
            'fc': self.fcd(),
            'eps_c': self.eps_c1,
            'eps_cu': self.eps_cu1,
            'Ec': self.Eci,
        }
