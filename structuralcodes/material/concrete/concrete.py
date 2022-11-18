"""Core implementation of the concrete material"""
import abc
import typing as t

# from structuralcodes.code import _CODE
from structuralcodes.code import _use_design_code

# To be done: rafactor in multiple files
from structuralcodes.code.mc2010 import mc2010
from structuralcodes.core.base import Material


class Concrete(Material):
    """The abstract concrete material."""

    _fck: float
    _existing: bool

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400,
        existing: t.Optional[bool] = False,
    ) -> None:
        '''Initializes an abstract concrete material'''
        name = name if name is not None else "Concrete"
        super().__init__(density=density, name=name)

        self._fck = abs(fck)
        self._existing = existing

    @property
    def fck(self) -> float:
        '''Returns fck in MPa'''
        return self._fck

    @fck.setter
    def fck(self, fck: float) -> None:
        '''Setter for fck (in MPa)'''
        self._fck = abs(fck)
        self._reset_attributes()

    @abc.abstractmethod
    def _reset_attributes(self):
        '''Each concrete should define its own _reset_attributes method
        This is because fck setting, reset the object arguments'''


def create_concrete(
    fck: float,
    name: t.Optional[str] = None,
    density: float = 2400.0,
    existing: bool = False,
    design_code: t.Optional[str] = None,
) -> t.Optional[Concrete]:
    """Create a concrete specifying the code"""
    # Get the code
    _code = _use_design_code(design_code)
    # Check if the code is a proper concrete code
    code = _code if 'concrete' in _code.__materials__ else None
    if code is None:
        raise Exception(
            'The design code is not set, either use '
            'structuralcodes.code.set_designcode, or provide a valid '
            'string in the function.'
        )
    # Create the proper concrete object
    if code.__name__ == 'structuralcodes.code.mc2010':
        return ConcreteMC2010(fck, name, density, existing)
    return None


# To be refactored in multiple files


class ConcreteMC2010(Concrete):
    '''Concrete implementation for MC 2010'''

    _fcm: t.Optional[float] = None
    _fctm: t.Optional[float] = None
    _fctkmin: t.Optional[float] = None
    _fctkmax: t.Optional[float] = None
    _Gf: t.Optional[float] = None
    _fcd: t.Optional[float] = None

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400.0,
        existing: bool = False,
    ):
        '''
        Initializes a new instance of Concrete for MC 2010

        :param float fck: Characteristic strength in MPa if concrete is not
            existing. Otherwise it is interpretedas fcm
        :param Optional[str] name: A descriptive name for concrete

        '''

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
        self._fcd = None

    def update_attributes(self, updated_attributes: dict) -> None:
        '''
        Function for updating the attributes specified in the input dictionary

        :param dict updated_attributes: the dictionary of parameters to be
            updated (not found parameters are skipped with a warning)
        '''
        for key, value in updated_attributes.items():
            if not hasattr(self, '_' + key):
                print(
                    f'WARNING: attribute {key} not found. Ignoring this entry'
                )
                str_list_keys = ''
                for k in updated_attributes.keys():
                    str_list_keys += k + ', '
                print(' Available keys: ' + str_list_keys)
                continue
            setattr(self, '_' + key, value)

    @property
    def fcm(self) -> float:
        '''Returns fcm in MPa'''
        if self._fcm is not None:
            return self._fcm
        return mc2010.fcm(self._fck)

    @fcm.setter
    def fcm(self, value: float):
        '''
        Sets a user defined value for fcm

        :param float value: the values of fcm
        :raises ValueError: if value is larger than fck'''
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
        '''Returns fctm in MPa'''
        if self._fctm is not None:
            return self._fctm
        return mc2010.fctm(self._fck)

    @fctm.setter
    def fctm(self, value: float):
        '''
        Sets a user defined value for fctm

        :param float value: the values of fctm
        '''
        # We could check and raise valueErrors if the ranges are "strange"
        # with respect to fck (e.g. if it is higher than x*fck)
        self._fctm = abs(value)

    @property
    def fctkmin(self) -> float:
        '''Returns fctkmin in MPa'''
        if self._fctkmin is not None:
            return self._fctkmin
        # TO DO: if fctm is not None: fctkmin is comput from that!
        return mc2010.fctkmin(self._fck)

    @fctkmin.setter
    def fctkmin(self, value: float):
        '''
        Sets a user defined value for fctkmin

        :param float value: the values of fctkmin
        '''
        # We could check and raise valueErrors if the ranges are "strange"
        # with respect to fck (e.g. if it is higher than x*fck)
        # or that fctkmax is > than fctm > fctkmin
        self._fctkmin = abs(value)

    @property
    def fctkmax(self) -> float:
        '''Returns fctkmax in MPa'''
        if self._fctkmax is not None:
            return self._fctkmax
        # TO DO: if fctm is not None: fctkmin is comput from that!
        return mc2010.fctkmax(self._fck)

    @fctkmax.setter
    def fctkmax(self, value: float):
        '''
        Sets a user defined value for fctkmax

        :param float value: the values of fctkmax
        '''
        # We could check and raise valueErrors if the ranges are "strange"
        # with respect to fck (e.g. if it is higher than x*fck)
        # or that fctkmax is > than fctm > fctkmin
        self._fctkmax = abs(value)

    @property
    def Gf(self) -> float:
        '''Returns Gf in N/m'''
        if self._Gf is not None:
            return self._Gf
        return mc2010.Gf(self._fck)

    @Gf.setter
    def Gf(self, value: float):
        '''
        Sets a user defined value for Gf

        :param float value: the values of Gf
        '''
        # We could check and raise valueErrors if the ranges are "strange"
        # with respect to fck (e.g. if it is higher than x*fck)
        # or that fctkmax is > than fctm > fctkmin
        self._Gf = abs(value)

    @property
    def fcd(self) -> float:
        '''Returns fcd in MPa'''
        if self._fcd is not None:
            return self._fcd
        return mc2010.fcd(self._fck)

    @fcd.setter
    def fcd(self, value: float):
        '''
        Sets a user defined value for fcd

        :param float value: the values of fcd
        '''
        # We could check and raise valueErrors if the values are "strange"
        self._fcd = abs(value)
