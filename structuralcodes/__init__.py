"""A Python package that contains models from structural design codes."""

import warnings

from . import codes, core, geometry, materials, sections
from .codes import get_design_codes, set_design_code, set_national_annex
from .core.errors import StructuralCodesWarning

__version__ = '0.6.0'

__all__ = [
    'set_design_code',
    'get_design_codes',
    'set_national_annex',
    'codes',
    'core',
    'materials',
    'geometry',
    'sections',
]

warnings.filterwarnings(action='always', category=StructuralCodesWarning)
