"""Core implementation of the concrete material"""
import typing as t

# from structuralcodes.code import _CODE
from structuralcodes.code import _use_design_code
from structuralcodes.core.base import DesignCode, Material


REQUIRED_FUNCTIONS = (
    'fcm',
    'fctm',
    'fctkmin',
    'fctkmax',
    'Gf',
    'fcd',
)

# TO BE DONE HERE: import Concrete Specific classes (MC, EN, etc.)?
# Create a concrete facotry to manage that?


class Concrete(Material):
    """Base class for all concete materials."""

    _code: DesignCode
    _fck: float
    _existing: bool
    _fcm: float
    _fctm: float
    _fctkmin: float
    _fctkmax: float
    _Gf: float
    _fcd: float

    def __init__(
        self, fck: float, design_code: t.Optional[str] = None
    ) -> None:
        super().__init__(name=f'C{fck}', density=2400)

        _code = _use_design_code(design_code)
        code = _code if 'concrete' in _code.materials else None

        if code is None:
            raise Exception(
                'The design code is not set, either use '
                'structuralcodes.code.set_designcode, or provide a valid '
                'string in the material constructor.'
            )

        self._code = code
        # fck is a positive number, take absolute value
        self._fck = abs(float(fck))
        # By the default concrete is a new concrete
        self._existing = False

        self._reset_attributes()

    def _reset_attributes(self):
        # Set all attributes to None, but anyhow check if the
        # function is defined in the design code
        for fun in REQUIRED_FUNCTIONS:
            if not hasattr(self._code, fun):
                raise Exception(
                    f'The function {fun} has not be implemented in code'
                )
            setattr(self, '_' + fun, None)

    def update_attributes(self, attrs_to_update: dict):
        '''Function for updating the attributes specified in the
        input dictionary.'''
        for key, value in attrs_to_update.items():
            if not hasattr(self, '_' + key):
                print(
                    f'WARNING: attribute {key} not found. Ignoring this entry'
                )
                continue
            setattr(self, '_' + key, value)

    @property
    def fck(self) -> float:
        '''Getter for property fck'''
        return self._fck

    @fck.setter
    def fck(self, value: float) -> float:
        '''Setting fck reset all attributes to the selected code'''
        self._fck = abs(value)

        self._reset_attributes()

    # All other properties call the proper function from the code
    @property
    def fcm(self):
        '''Getter for property fcm'''
        if self._fcm is not None:
            return self._fcm
        else:
            return getattr(self._code, 'fcm')(self._fck)

    @fcm.setter
    def fcm(self, value):
        '''Setting fcm: this substitutes the value computed by the current
        code.
        To reset, set property fck again'''
        self._fcm = abs(value)
        # Probably here we could set _fck based on _fcm (useful for
        #  existing structures)

    @property
    def fctm(self):
        '''Getter for property fctm'''
        if self._fctm is not None:
            return self._fctm
        else:
            return getattr(self._code, 'fctm')(self._fck)

    @fctm.setter
    def fctm(self, value):
        '''Setting fctm: this substitutes the value computed by the
        current code.
        To reset, set property fck again'''
        self._fctm = value

    @property
    def fctkmin(self):
        '''Getter for property fctkmin'''
        if self._fctkmin is not None:
            return self._fctkmin
        else:
            return getattr(self._code, 'fctkmin')(self._fck)

    @fctkmin.setter
    def fctkmin(self, value):
        '''Setting fctkmin: this substitutes the value computed by the
         current code.
        To reset, set property fck again'''
        # Here we could do some checks on value before assigning it.
        # For instance check fctkmin that has legit values
        self._fctkmin = value

    @property
    def fctkmax(self):
        '''Getter for property fctkmax'''
        if self._fctkmax is not None:
            return self._fctkmax
        else:
            return getattr(self._code, 'fctkmax')(self._fck)

    @fctkmax.setter
    def fctkmax(self, value):
        '''Setting fctkmin: this substitutes the value computed by the
        current code.
        To reset, set property fck again'''
        # Here we could do some checks on value before assigning it.
        # For instance check fctkmax
        self._fctkmax = value

    @property
    def fcd(self):
        '''Getter for property fcd'''
        if self._fcd is not None:
            return self._fcd
        else:
            return getattr(self._code, 'fcd')(self._fck)

    @fcd.setter
    def fcd(self, value):
        '''Setting fctkmin: this substitutes the value computed by the
        current code.
        To reset, set property fck again'''
        # Here we could do some checks on value before assigning it.
        # For instance check
        self._fcd = value

    @property
    def existing(self):
        '''Getter for property existing'''
        return self._existing

    @existing.setter
    def existing(self, value: bool):
        '''Description here'''
        FC = 1.0  # For now it is here, but we should think if
        # putting it in the code, or in concrete...
        self._existing = value
        self._fcd = getattr(self._code, 'fcd')(
            self._fck, existing=value, FC=FC
        )
