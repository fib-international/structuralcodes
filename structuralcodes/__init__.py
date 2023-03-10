"""A Python package that contains models from structural design codes"""
from .codes import set_design_code, get_design_codes, set_national_annex

from . import materials
from . import core
from . import codes

__version__ = ''

__all__ = [
    'set_design_code',
    'get_design_codes',
    'set_national_annex',
    'codes',
    'core',
    'materials',
]
