"""A Python package that contains models from structural design codes"""
from .code import set_design_code, set_national_annex
from .code.mc2010.mc2010 import MC2010
from .core.concrete import Concrete

__version__ = ''

__all__ = [
    'set_design_code',
    'set_national_annex',
    'Concrete',
    'MC2010',
]
