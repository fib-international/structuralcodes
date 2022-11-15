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


class Concrete(Material):
    """The concrete material."""

    _code: DesignCode
    _fck: float
    _fcm: float
    _fctm: float
    _fctkmin: float
    _fctkmax: float
    _Gf: float
    _fcd: float
    _existing: bool

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
        self._fck = abs(float(fck))
        self._existing = False

        # Set attributes from the design code
        for fun in REQUIRED_FUNCTIONS:
            if not hasattr(code, fun):
                raise Exception(
                    f'The function {fun} has not be implemented in code'
                )
            setattr(self, '_' + fun, getattr(code, fun)(fck))
            # Note: according to this approach everything is function of fck
            #  -> possible weak point?

    @property
    def fck(self) -> float:
        '''Getter for property fck'''
        return self._fck

    @fck.setter
    def fck(self, value: float) -> float:
        '''Setting fck reset all attributes to the selected code'''
        # check on fck value to be a legit one -> this depends on the standard?
        if value < 0:
            value *= -1
            # alternatively raise ValueError("Raise exceptions on values?")
        self._fck = value

        # Set attributes from the design code
        for fun in REQUIRED_FUNCTIONS:
            # This check maybe is not needed anymore?
            if not hasattr(self._code, fun):
                raise Exception(
                    f'The function {fun} has not be implemented in code'
                )
            setattr(self, '_' + fun, getattr(self._code, fun)(self._fck))
            # Note: according to this approach everything is function of fck
            #  -> possible weak point?

    @property
    def fcm(self):
        '''Getter for property fcm'''
        return self._fcm

    @fcm.setter
    def fcm(self, value):
        '''Setting fcm: this substitutes the value computed by the current
        code.
        To reset, set property fck again'''
        # Here we could do some checks on value before assigning it.
        # For instance check
        self._fcm = value
        # Probably here we could set _fck based on _fcm (useful for
        #  existing structures)

    @property
    def fctm(self):
        '''Getter for property fctm'''
        return self._fctm

    @fctm.setter
    def fctm(self, value):
        '''Setting fctm: this substitutes the value computed by the
        current code.
        To reset, set property fck again'''
        # Here we could do some checks on value before assigning it.
        # For instance check
        self._fctm = value

    @property
    def fctkmin(self):
        '''Getter for property fctkmin'''
        return self._fctkmin

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
        return self._fctkmax

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
        return self._fcd

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
        self.fcd = getattr(self._code, 'fcd')(self._fck, existing=value, FC=FC)
