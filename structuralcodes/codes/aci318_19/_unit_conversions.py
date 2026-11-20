"""US customary and SI unit conversions for ACI 318-19 helpers."""

PSI_TO_MPA = 0.006894757293168361
KSI_TO_MPA = 6.894757293168361
INCH_TO_MM = 25.4
LB_FORCE_TO_N = 4.4482216152605
KIP_TO_N = 1000.0 * LB_FORCE_TO_N
LB_PER_CUBIC_FOOT_TO_KG_PER_CUBIC_METER = 16.01846337396014


def psi_to_mpa(psi: float) -> float:
    """Convert stress from psi to MPa."""
    return psi * PSI_TO_MPA


def mpa_to_psi(mpa: float) -> float:
    """Convert stress from MPa to psi."""
    return mpa / PSI_TO_MPA


def ksi_to_mpa(ksi: float) -> float:
    """Convert stress from ksi to MPa."""
    return ksi * KSI_TO_MPA


def mpa_to_ksi(mpa: float) -> float:
    """Convert stress from MPa to ksi."""
    return mpa / KSI_TO_MPA


def pcf_to_kg_per_m3(pcf: float) -> float:
    """Convert density from lb/ft3 to kg/m3."""
    return pcf * LB_PER_CUBIC_FOOT_TO_KG_PER_CUBIC_METER


def kg_per_m3_to_pcf(kg_per_m3: float) -> float:
    """Convert density from kg/m3 to lb/ft3."""
    return kg_per_m3 / LB_PER_CUBIC_FOOT_TO_KG_PER_CUBIC_METER


def in_to_mm(inches: float) -> float:
    """Convert length from inches to mm."""
    return inches * INCH_TO_MM


def mm_to_in(mm: float) -> float:
    """Convert length from mm to inches."""
    return mm / INCH_TO_MM


def in2_to_mm2(square_inches: float) -> float:
    """Convert area from in2 to mm2."""
    return square_inches * INCH_TO_MM**2


def mm2_to_in2(square_mm: float) -> float:
    """Convert area from mm2 to in2."""
    return square_mm / INCH_TO_MM**2


def in3_to_mm3(cubic_inches: float) -> float:
    """Convert volume from in3 to mm3."""
    return cubic_inches * INCH_TO_MM**3


def mm3_to_in3(cubic_mm: float) -> float:
    """Convert volume from mm3 to in3."""
    return cubic_mm / INCH_TO_MM**3


def in4_to_mm4(inches_fourth: float) -> float:
    """Convert second moment of area from in4 to mm4."""
    return inches_fourth * INCH_TO_MM**4


def mm4_to_in4(mm_fourth: float) -> float:
    """Convert second moment of area from mm4 to in4."""
    return mm_fourth / INCH_TO_MM**4


def lb_to_n(pounds: float) -> float:
    """Convert force from lbf to N."""
    return pounds * LB_FORCE_TO_N


def n_to_lb(newtons: float) -> float:
    """Convert force from N to lbf."""
    return newtons / LB_FORCE_TO_N


def kip_to_n(kips: float) -> float:
    """Convert force from kip to N."""
    return kips * KIP_TO_N


def n_to_kip(newtons: float) -> float:
    """Convert force from N to kip."""
    return newtons / KIP_TO_N


def lb_in_to_nmm(lb_in: float) -> float:
    """Convert moment from lb-in to Nmm."""
    return lb_in * LB_FORCE_TO_N * INCH_TO_MM


def nmm_to_lb_in(nmm: float) -> float:
    """Convert moment from Nmm to lb-in."""
    return nmm / (LB_FORCE_TO_N * INCH_TO_MM)


def kip_in_to_nmm(kip_in: float) -> float:
    """Convert moment from kip-in to Nmm."""
    return kip_in * KIP_TO_N * INCH_TO_MM


def nmm_to_kip_in(nmm: float) -> float:
    """Convert moment from Nmm to kip-in."""
    return nmm / (KIP_TO_N * INCH_TO_MM)
