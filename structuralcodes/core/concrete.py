"""Core implementation of the concrete material"""
import typing as t

# from structuralcodes.code import _CODE
from structuralcodes.code import _use_design_code
from structuralcodes.core.base import ConcreteDesignCode, Material, DesignCode


REQUIRED_FUNCTIONS = (
    'fcm',
    'fctm',
    'fctkmin',
    'fctkmax',
    'Gf',
)


class Concrete(Material):
    """The concrete material."""
    _code: ConcreteDesignCode
    _fck: float
    _fcm: float
    _fctm: float
    _fctkmin: float
    _fctkmax: float
    _Gf: float

    def __init__(self, fck: float, design_code: t.Optional[str] = None) -> None:
        super().__init__()

        _code = _use_design_code(design_code)
        code = _code if isinstance(_code, DesignCode) else None

        if code is None:
            raise Exception(
                'The design code is not set, either use '
                'structuralcodes.code.set_designcode, or provide a valid '
                'string in the material constructor.'
            )

        # check on fck value to be a legit one -> this depends on the standard?
        if fck < 0:
            fck *= -1

        self._code = code
        self._fck = fck

        # Set attributes from the design code
        for fun in REQUIRED_FUNCTIONS:
            if not hasattr(code,fun):
                raise Exception(
                    'The function {} has not be implemented in code'.format(
                        fun
                    )
                )
            setattr(self, '_'+fun, getattr(code,fun)(fck))
            # Note: according to this approach everything is function of fck -> possible weak point

    @property
    def fck(self):
        return self._fck

    @fck.setter
    def fck(self, value):
        ''' Setting fck reset all attributes to the selected code'''
        # check on fck value to be a legit one -> this depends on the standard?
        if value <0:
            fck *= -1
            # raise ValueError("Raise exceptions on values?")
        self._fck = value

        # Set attributes from the design code
        for fun in REQUIRED_FUNCTIONS:
            # This check maybe is not needed anymore?
            if not hasattr(self._code,fun):
                raise Exception(
                    'The function {} has not be implemented in code'.format(
                        fun.__name__
                    )
                )
            setattr(self, '_'+fun, self._code.fun(self._fck))
            # Note: according to this approach everything is function of fck -> possible weak point
    
    @property
    def fcm(self):
        return(self._fcm)
    
    @fcm.setter
    def fcm(self, value):
        '''Setting fcm: this substitutes the value computed by the current code.
        To reset, set property fck again'''
        # Here we could do some checks on value before assigning it. For instance check 
        self._fcm = value
        #Probably here we could set _fck based on _fcm (useful for existing structures)
    
    @property 
    def fctm(self):
        return(self._fctm)
    
    @fctm.setter
    def fctm(self,value):
        '''Setting fctm: this substitutes the value computed by the current code.
        To reset, set property fck again'''
        # Here we could do some checks on value before assigning it. For instance check 
        self._fctm = value

    @property 
    def fctkmin(self):
        return(self._fctkmin)
    
    @fctkmin.setter
    def fctkmin(self,value):
        '''Setting fctkmin: this substitutes the value computed by the current code.
        To reset, set property fck again'''
        # Here we could do some checks on value before assigning it. For instance check 
        self._fctkmin = value

    @property 
    def fctkmax(self):
        return(self._fctkmax)
    
    @fctkmax.setter
    def fctkmax(self,value):
        '''Setting fctkmin: this substitutes the value computed by the current code.
        To reset, set property fck again'''
        # Here we could do some checks on value before assigning it. For instance check 
        self._fctkmax = value
