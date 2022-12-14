"""A Python package that contains models from structural design codes"""
from .code import set_design_code, get_design_codes, set_national_annex

from .code import mc2010

__version__ = ''

__all__ = [
    'set_design_code',
    'get_design_codes',
    'set_national_annex',
    'mc2010',
]
