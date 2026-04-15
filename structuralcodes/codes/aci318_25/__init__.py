"""ACI 318-25: Building Code Requirements for Structural Concrete."""

import typing as t

from ._units import (
    FT_TO_MM,
    IN_TO_MM,
    KGM3_TO_PCF,
    KIP_TO_N,
    KSI_TO_MPA,
    LBF_TO_N,
    MM_TO_FT,
    MM_TO_IN,
    MPA_TO_KSI,
    MPA_TO_PSI,
    N_TO_KIP,
    N_TO_LBF,
    PCF_TO_KGM3,
    PSF_TO_KPA,
    PSF_TO_PA,
    PSI_TO_MPA,
)

__all__ = [
    'PSI_TO_MPA',
    'KSI_TO_MPA',
    'MPA_TO_PSI',
    'MPA_TO_KSI',
    'IN_TO_MM',
    'FT_TO_MM',
    'MM_TO_IN',
    'MM_TO_FT',
    'LBF_TO_N',
    'KIP_TO_N',
    'N_TO_LBF',
    'N_TO_KIP',
    'PSF_TO_PA',
    'PSF_TO_KPA',
    'PCF_TO_KGM3',
    'KGM3_TO_PCF',
]

__title__: str = 'ACI 318-25'
__year__: str = '2025'
__materials__: t.Tuple[str, ...] = ('concrete', 'reinforcement')
