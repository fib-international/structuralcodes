"""The concrete class for EC2 2004 Concrete Material."""

import typing as t
import warnings

from structuralcodes.codes import ec2_2004

from ._concrete import Concrete


class ConcreteEC2_2004(Concrete):  # noqa: N801
    """Concrete implementation for EC2 2004."""

    _fcm: t.Optional[float] = None
    _fctm: t.Optional[float] = None
    _fctk_5: t.Optional[float] = None
    _fctk_95: t.Optional[float] = None
    _Ecm: t.Optional[float] = None

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400,
        **kwargs,
    ) -> None:
        """Initializes a new instance of Concrete for EC2 2004.

        Args:
            fck (float): Characteristic strength in MPa if concrete is not
                existing.

        Keyword Args:
            name (str): A descriptive name for concrete
            density (float): Density of material in kg/m3 (default: 2400)
            existing (bool): The material is of an existing structure
                (default: False)
        """
        del kwargs
        if name is None:
            name = f'C{round(fck):d}'
        super().__init__(fck=fck, name=name, density=density, existing=False)

    def _reset_attributes(self):
        self._fcm = None
        self._fctm = None
        self._fctk_5 = None
        self._fctk_95 = None
        self._Ecm = None

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
        if self._fctm is None:
            self._fctm = ec2_2004.fctm(self._fck)
        return self._fctm

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
    def fctk_5(self) -> float:
        """Returns fctk_5 in MPa.

        Returns:
            float: The lower bound tensile strength in MPa
        """
        if self._fctk_5 is not None:
            return self._fctk_5
        return ec2_2004.fctk_5(self.fctm)

    @fctk_5.setter
    def fctk_5(self, value: float):
        """Sets a user defined value for fctk_5.

        Args:
            value (float): the value of fctk_5 in MPa
        """
        self._fctk_5 = abs(value)

    @property
    def fctk_95(self) -> float:
        """Returns fctk_95 in MPa.

        Returns:
            float: The upper bound tensile strength in MPa
        """
        if self._fctk_95 is not None:
            return self._fctk_95

        return ec2_2004.fctk_95(self.fctm)

    @fctk_95.setter
    def fctk_95(self, value: float):
        """Sets a user defined value for fctk_95.

        Args:
            value (float): the value of fctk_95 in MPa
        """
        self._fctk_95 = abs(value)

    @property
    def Ecm(self) -> float:
        """Returns Ecm in MPa.

        Returns:
            float: The upper bound tensile strength in MPa
        """
        if self._Ecm is not None:
            return self._Ecm

        return ec2_2004.Ecm(self.fcm)

    @Ecm.setter
    def Ecm(self, value: float):
        """Sets a user defined value for Ecm.

        Args:
            value (float): the value of Ecm in MPa
        """
        self._Ecm = abs(value)
