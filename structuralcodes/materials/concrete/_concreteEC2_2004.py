"""The concrete class for EC2 2004 Concrete Material."""

import typing as t
import warnings

from structuralcodes.codes import ec2_2004

from ._concrete import Concrete


class ConcreteEC2_2004(Concrete):  # noqa: N801
    """Concrete implementation for EC2 2004."""

    # computed attributes
    _fcm: t.Optional[float] = None
    _fctm: t.Optional[float] = None
    _fctk_5: t.Optional[float] = None
    _fctk_95: t.Optional[float] = None
    _Ecm: t.Optional[float] = None
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
        density: float = 2400,
        gamma_c: t.Optional[float] = None,
        alpha_cc: t.Optional[float] = None,
        **kwargs,
    ) -> None:
        """Initializes a new instance of Concrete for EC2 2004.

        Arguments:
            fck (float): Characteristic strength in MPa if concrete is not
                existing.

        Keyword Arguments:
            name (str): A descriptive name for concrete.
            density (float): Density of material in kg/m3 (default: 2400).
            gamma_c (float, optional): partial factor of concrete (default is
                1.5).
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
            existing=False,
            gamma_c=gamma_c,
        )
        self._alpha_cc = alpha_cc

    def _reset_attributes(self):
        self._fcm = None
        self._fctm = None
        self._fctk_5 = None
        self._fctk_95 = None
        self._Ecm = None
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
        if self._fcm is None:
            self._fcm = ec2_2004.fcm(self._fck)
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
    def fctm(self) -> float:
        """Returns fctm in MPa.

        Returns:
            float: The mean tensile strength in MPa.
        """
        if self._fctm is None:
            self._fctm = ec2_2004.fctm(self._fck)
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
    def fctk_5(self) -> float:
        """Returns fctk_5 in MPa.

        Returns:
            float: The lower bound tensile strength in MPa.
        """
        if self._fctk_5 is not None:
            return self._fctk_5
        return ec2_2004.fctk_5(self.fctm)

    @fctk_5.setter
    def fctk_5(self, value: float):
        """Sets a user defined value for fctk_5.

        Arguments:
            value (float): The value of fctk_5 in MPa.
        """
        self._fctk_5 = abs(value)

    @property
    def fctk_95(self) -> float:
        """Returns fctk_95 in MPa.

        Returns:
            float: The upper bound tensile strength in MPa.
        """
        if self._fctk_95 is not None:
            return self._fctk_95

        return ec2_2004.fctk_95(self.fctm)

    @fctk_95.setter
    def fctk_95(self, value: float):
        """Sets a user defined value for fctk_95.

        Arguments:
            value (float): The value of fctk_95 in MPa.
        """
        self._fctk_95 = abs(value)

    @property
    def Ecm(self) -> float:
        """Returns Ecm in MPa.

        Returns:
            float: The upper bound tensile strength in MPa.
        """
        if self._Ecm is not None:
            return self._Ecm

        return ec2_2004.Ecm(self.fcm)

    @Ecm.setter
    def Ecm(self, value: float):
        """Sets a user defined value for Ecm.

        Arguments:
            value (float): The value of Ecm in MPa.
        """
        self._Ecm = abs(value)

    def fcd(self) -> float:
        """Return the design compressive strength in MPa.

        Returns:
            float: The design compressive strength of concrete in MPa.
        """
        # This method should perhaps become a property, but is left as a method
        # for now, to be consistent with other concretes.
        return ec2_2004.fcd(
            self.fck, alpha_cc=self.alpha_cc, gamma_c=self.gamma_c
        )

    @property
    def gamma_c(self) -> float:
        """The partial factor for concrete."""
        # Here we should implement the interaction with the globally set
        # national annex. For now, we simply return the default value.
        return self._gamma_c or 1.5

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
        self._eps_c1 = self._eps_c1 or ec2_2004.eps_c1(self.fcm)
        return self._eps_c1

    @eps_c1.setter
    def eps_c1(self, value: float):
        """Sets a user defined value for strain at peak strength for Sargin
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
        self._eps_cu1 = self._eps_cu1 or ec2_2004.eps_cu1(self.fcm)
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
        self._k_sargin = self._k_sargin or ec2_2004.k_sargin(
            Ecm=self.Ecm,
            fcm=self.fcm,
            eps_c1=self.eps_c1,
        )
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
        self._eps_c2 = self._eps_c2 or ec2_2004.eps_c2(self.fck)
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
        self._eps_cu2 = self._eps_cu2 or ec2_2004.eps_cu2(self.fck)
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
            float: The exponent for Parabola-rectangle law.
        """
        self._n_parabolic_rectangular = (
            self._n_parabolic_rectangular
            or ec2_2004.n_parabolic_rectangular(self.fck)
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
        self._eps_c3 = self._eps_c3 or ec2_2004.eps_c3(self.fck)
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
        self._eps_cu3 = self._eps_cu3 or ec2_2004.eps_cu3(self.fck)
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
        return {'E': self.Ecm}

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
            'Ec': self.Ecm,
        }
