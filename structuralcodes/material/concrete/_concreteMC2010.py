"""The concrete class for Model Code 2020 Concrete Material"""
import typing as t
import warnings

from structuralcodes.codes import mc2010
from ._concrete import Concrete


class ConcreteMC2010(Concrete):
    """Concrete implementation for MC 2010"""

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
        existing: bool = False,
    ):
        """Initializes a new instance of Concrete for MC 2010

        Args:
            fck (float): Characteristic strength in MPa if concrete is not
                existing.

        Keyword Args:
            name (str): A descriptive name for concrete
            density (float): Density of material in kg/m3 (default: 2400)
            existing (bool): The material is of an existing structure
                (default: False)
        """

        if name is None:
            name = f'C{round(fck):d}'
        super().__init__(
            fck=fck, name=name, density=density, existing=existing
        )

    def _reset_attributes(self):
        self._fcm = None
        self._fctm = None
        self._fctkmin = None
        self._fctkmax = None
        self._Gf = None

    def update_attributes(self, updated_attributes: dict) -> None:
        """Function for updating the attributes specified in the input
        dictionary

        Args:
            updated_attributes (dict): the dictionary of parameters to be
                updated (not found parameters are skipped with a warning)
        """
        for key, value in updated_attributes.items():
            if not hasattr(self, '_' + key):
                str_list_keys = ''
                for k in updated_attributes.keys():
                    str_list_keys += k + ', '
                str_warn = (
                    f'WARNING: attribute {key} not found. Ignoring the entry.'
                )
                str_warn += '\nAvailable keys: ' + str_list_keys
                warnings.warn(str_warn)
                continue
            setattr(self, '_' + key, value)

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
        """Sets a user defined value for fcm

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
        """Returns fctm in MPa

        Returns:
            float: The mean tensile strength in MPa
        """
        if self._fctm is not None:
            return self._fctm
        return mc2010.fctm(self._fck)

    @fctm.setter
    def fctm(self, value: float):
        """Sets a user defined value for fctm

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
        """Returns fctkmin in MPa

        Returns:
            float: The lower bound tensile strength in MPa
        """
        if self._fctkmin is not None:
            return self._fctkmin

        return mc2010.fctkmin(self.fctm)

    @fctkmin.setter
    def fctkmin(self, value: float):
        """Sets a user defined value for fctkmin

        Args:
            value (float): the value of fctkmin in MPa
        """
        self._fctkmin = abs(value)

    @property
    def fctkmax(self) -> float:
        """Returns fctkmax in MPa

        Returns:
            float: The upper bound tensile strength in MPa
        """
        if self._fctkmax is not None:
            return self._fctkmax

        return mc2010.fctkmax(self.fctm)

    @fctkmax.setter
    def fctkmax(self, value: float):
        """Sets a user defined value for fctkmax

        Args:
            value (float): the value of fctkmax in MPa
        """
        self._fctkmax = abs(value)

    @property
    def Gf(self) -> float:
        """Fracture energy of concrete

        Returns:
            float: The fracture energy in N/m
        """
        if self._Gf is not None:
            return self._Gf
        return mc2010.Gf(self._fck)

    @Gf.setter
    def Gf(self, value: float):
        """Sets a user defined value for fracture energy Gf

        Args:
            value (float): the value of Gf in N/m
        """
        self._Gf = abs(value)
