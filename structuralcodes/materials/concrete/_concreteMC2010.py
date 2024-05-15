"""The concrete class for Model Code 2020 Concrete Material."""

import typing as t
import warnings

from structuralcodes.codes import mc2010

from ._concrete import Concrete


class ConcreteMC2010(Concrete):
    """Concrete implementation for MC 2010."""

    _fcm: t.Optional[float] = None
    _fctm: t.Optional[float] = None
    _fctkmin: t.Optional[float] = None
    _fctkmax: t.Optional[float] = None
    _Gf: t.Optional[float] = None

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400.0,
        gamma_c: t.Optional[float] = None,
        existing: bool = False,
    ):
        """Initializes a new instance of Concrete for MC 2010.

        Args:
            fck (float): Characteristic strength in MPa if concrete is not
                existing.

        Keyword Args:
            name (Optional(str)): A descriptive name for concrete.
            density (float): Density of material in kg/m3 (default: 2400).
            gamma_c (Optional(float)): The partial factor for concrete.
            existing (bool): The material is of an existing structure
                (default: False).
        """
        if name is None:
            name = f'C{round(fck):d}'
        super().__init__(
            fck=fck,
            name=name,
            density=density,
            existing=existing,
            gamma_c=gamma_c,
        )

    def _reset_attributes(self):
        self._fcm = None
        self._fctm = None
        self._fctkmin = None
        self._fctkmax = None
        self._Gf = None

    @property
    def fcm(self) -> float:
        """Returns fcm in MPa.

        Returns:
            float: The mean compressive strength in MPa.
        """
        if self._fcm is not None:
            return self._fcm
        return mc2010.fcm(self._fck)

    @fcm.setter
    def fcm(self, value: float):
        """Sets a user defined value for fcm.

        Args:
            value (float): the value of fcm in MPa

        Raises:
            ValueError: if value is lower than fck
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
            float: The mean tensile strength in MPa
        """
        if self._fctm is not None:
            return self._fctm
        return mc2010.fctm(self._fck)

    @fctm.setter
    def fctm(self, value: float):
        """Sets a user defined value for fctm.

        Args:
            value (float): the value of fctm in MPa
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
            float: The lower bound tensile strength in MPa
        """
        if self._fctkmin is not None:
            return self._fctkmin

        return mc2010.fctkmin(self.fctm)

    @fctkmin.setter
    def fctkmin(self, value: float):
        """Sets a user defined value for fctkmin.

        Args:
            value (float): the value of fctkmin in MPa
        """
        self._fctkmin = abs(value)

    @property
    def fctkmax(self) -> float:
        """Returns fctkmax in MPa.

        Returns:
            float: The upper bound tensile strength in MPa
        """
        if self._fctkmax is not None:
            return self._fctkmax

        return mc2010.fctkmax(self.fctm)

    @fctkmax.setter
    def fctkmax(self, value: float):
        """Sets a user defined value for fctkmax.

        Args:
            value (float): the value of fctkmax in MPa
        """
        self._fctkmax = abs(value)

    @property
    def Gf(self) -> float:
        """Fracture energy of concrete.

        Returns:
            float: The fracture energy in N/m
        """
        if self._Gf is not None:
            return self._Gf
        return mc2010.Gf(self._fck)

    @Gf.setter
    def Gf(self, value: float):
        """Sets a user defined value for fracture energy Gf.

        Args:
            value (float): the value of Gf in N/m
        """
        self._Gf = abs(value)

    @property
    def gamma_c(self) -> float:
        """The partial factor for concrete."""
        return self._gamma_c or 1.5

    def fcd(self, alpha_cc: float) -> float:
        """Calculate the design compressive strength.

        Args:
            alpha_cc (float): A factor for considering long-term effects on the
                strength, and effects that arise from the way the load is
                applied.

        Returns:
            float: The design compressive strength of concrete in MPa
        """
        return mc2010.fcd(self.fck, alpha_cc=alpha_cc, gamma_c=self.gamma_c)
