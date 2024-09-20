"""A Python package that contains models from structural design codes."""

from . import codes, core, geometry, materials, sections
from .codes import get_design_codes, set_design_code, set_national_annex

__version__ = '0.0.2'

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
