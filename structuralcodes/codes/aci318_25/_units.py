"""Unit conversion constants for ACI 318-25.

These constants convert between US customary units and SI units as used
in ACI 318-25: Building Code Requirements for Structural Concrete.
"""

# Pressure / stress conversions
PSI_TO_MPA: float = 0.00689476
KSI_TO_MPA: float = 6.89476
MPA_TO_PSI: float = 145.038
MPA_TO_KSI: float = 0.145038

# Length conversions
IN_TO_MM: float = 25.4
FT_TO_MM: float = 304.8
MM_TO_IN: float = 1 / 25.4
MM_TO_FT: float = 1 / 304.8

# Force conversions
LBF_TO_N: float = 4.44822
KIP_TO_N: float = 4448.22
N_TO_LBF: float = 1 / 4.44822
N_TO_KIP: float = 1 / 4448.22

# Distributed load / density conversions
PSF_TO_PA: float = 47.8803
PSF_TO_KPA: float = 0.0478803
PCF_TO_KGM3: float = 16.0185
KGM3_TO_PCF: float = 1 / 16.0185
