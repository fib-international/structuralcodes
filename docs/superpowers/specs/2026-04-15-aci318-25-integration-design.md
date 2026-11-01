# ACI 318-25 Integration into structuralcodes

**Date:** 2026-04-15
**Status:** Approved
**Scope:** One-way slab flexural and shear design (4000 psi concrete, Gr 60 rebar)

## Overview

Integrate ACI 318-25 (Building Code Requirements for Structural Concrete) into the
structuralcodes library. The design follows Approach 3 ("Thin Adapter, Thick Code Module"):
minimal material classes that plug into the existing geometry/section/integrator pipeline,
with the real ACI intelligence living in a self-contained code module.

Two design paths are supported:
- **Path A — Closed-form ACI equations:** Whitney stress block hand-calc functions for
  standard design (the typical engineering workflow).
- **Path B — Section-integrator-driven analysis:** Existing Marin/Fiber integrators with
  ACI-appropriate constitutive laws for non-standard situations (interaction diagrams,
  biaxial bending, moment-curvature).

### Key Architectural Decisions

1. **No changes to existing base classes, sections, geometry, or integrators.** The existing
   pipeline is consumed as-is.
2. **Safety philosophy:** ACI uses LRFD — material strengths are unreduced (`gamma_c=1.0`,
   `gamma_s=1.0`), and strength reduction factors (phi) are applied to member capacity.
   Phi logic is centralized in a dedicated module. Design functions return nominal strengths;
   phi is applied externally by the caller.
3. **Units:** SI internally (MPa, mm, N, N-mm) consistent with the library. A thin
   unit-conversion utility provides constants (`PSI_TO_MPA`, `IN_TO_MM`, etc.).
4. **Constitutive laws stay generic and shared.** Code-specific design idealizations
   (Whitney block, ACI strain limits) are parameterized through the material class dunder
   methods, following the existing factory pattern.
5. **Independent implementation.** Not dependent on PR #343 (gabe-kafka's ACI 318-19
   material layer). Can be reconciled later if both land upstream.

## 1. Code Module Structure

```
structuralcodes/codes/aci318_25/
├── __init__.py                            # Registry: __title__, __year__, __materials__, public API
├── _concrete_material_properties.py       # Ch. 19: Ec, fr, beta1, eps_cu, alpha1, fct, lambda_factor
├── _reinforcement_material_properties.py  # Ch. 20: Es, fy, grade lookup (ASTM A615 Gr 40/60/80/100)
├── _strength_reduction.py                 # Ch. 21, Table 21.2.1/21.2.2: phi factors, strain-based transition
├── _flexure.py                            # Ch. 22.2-22.3: Mn, As_required, As_min, As_max, c/d checks
├── _shear.py                              # Ch. 22.5: Vc (detailed + simplified), Vs, Vn, Av_min
├── _one_way_slab.py                       # Ch. 7: min thickness, shrinkage steel, bar spacing, critical section
└── _units.py                              # Unit conversion constants (PSI_TO_MPA, IN_TO_MM, etc.)
```

### Design principles

- All functions are pure — explicit parameters in, values out. No global state access.
- All units are SI (MPa, mm, mm2, N, N-mm) to match the library.
- Each function documents the ACI 318-25 section/equation it implements.
- Functions return nominal strengths. The phi module is called separately.

### Registration

In `codes/__init__.py`, add to imports and `_DESIGN_CODES`:
```python
_DESIGN_CODES = {
    'aci318_25': aci318_25,  # added
    'mc2010': mc2010,
    ...
}
```

Module metadata in `codes/aci318_25/__init__.py`:
```python
__title__: str = 'ACI 318-25'
__year__: str = '2025'
__materials__: tuple = ('concrete', 'reinforcement')
```

## 2. Material Classes

### ConcreteACI318_25

Inherits from `Concrete`. Satisfies the abstract interface while exposing ACI-native
properties.

**File:** `materials/concrete/_concreteACI318_25.py`

```python
class ConcreteACI318_25(Concrete):
    """ACI 318-25 concrete material.

    Uses LRFD philosophy — material strengths are unreduced. Safety is applied
    at member capacity level via strength reduction factors phi (Ch. 21).

    The gamma_c property returns 1.0 to satisfy the Concrete base class
    interface. ACI 318 does not use material partial factors. The fcd() method
    returns alpha1 * f'c (= 0.85 * f'c), which is the stress intensity used
    in the Whitney equivalent rectangular stress block, not a gamma-reduced
    design strength in the Eurocode sense.
    """
```

Properties:

| Property | Source | Notes |
|----------|--------|-------|
| `fck` | Input | Specified compressive strength f'c (MPa) |
| `fc` | Alias for `fck` | ACI notation convenience |
| `gamma_c` | `1.0` | Satisfies base class; ACI has no material partial factor |
| `fcd()` | `alpha1 * fc / gamma_c` | = 0.85 * f'c; feeds constitutive law factory |
| `Ec` | `aci318_25.Ec(fc, wc)` | Table 19.2.2.1 |
| `fr` | `aci318_25.fr(fc, lambda_s)` | Eq. 19.2.3.1 (modulus of rupture) |
| `fct` | `aci318_25.fct(fc, lambda_s)` | Sec. 19.2.4.3 (splitting tensile) |
| `beta1` | `aci318_25.beta1(fc)` | Table 22.2.2.4.3 |
| `alpha1` | `0.85` | Sec. 22.2.2.4.1 |
| `eps_cu` | `0.003` | Sec. 22.2.2.1 |
| `lambda_factor` | `aci318_25.lambda_factor(type)` | Table 19.2.4.2 |

Constitutive law dunder methods:

| Method | Returns |
|--------|---------|
| `__elastic__()` | `{'E': self.Ec}` |
| `__parabolarectangle__()` | `{'fc': self.fcd(), 'eps_0': 0.002, 'eps_u': 0.003, 'n': 2}` |
| `__bilinearcompression__()` | `{'fc': self.fcd(), 'eps_c': 0.002, 'eps_cu': 0.003}` |
| `__whitneyblock__()` | `{'fc': self.alpha1 * self.fc, 'beta1': self.beta1, 'eps_cu': 0.003}` |

Notes:
- `eps_0=0.002` for parabola-rectangle is the Hognestad peak strain.
- `eps_u=0.003` enforces the ACI ultimate strain assumption.
- `fcd()` returns `0.85 * f'c`, which feeds correctly into both Whitney and parabola-rectangle.

### ReinforcementACI318_25

**File:** `materials/reinforcement/_reinforcementACI318_25.py`

```python
class ReinforcementACI318_25(Reinforcement):
    """ACI 318-25 reinforcement material.

    Strengths are unreduced (gamma_s=1.0). ACI applies strength reduction
    factors at member capacity level, not material level.

    Supports ASTM A615 grades: 40, 60, 80, 100.
    """
```

Properties:

| Property | Value |
|----------|-------|
| `gamma_s` | `1.0` |
| `fyd()` | `fy` (unreduced) |
| `ftd()` | `fu` (unreduced) |
| `epsud()` | `eps_su` |

Convenience constructor:
```python
@classmethod
def from_grade(cls, grade='60') -> 'ReinforcementACI318_25':
    """Create from ASTM A615 grade. Looks up fy, fu from grade table."""
```

### Factory Registration

In `materials/concrete/__init__.py`, add `ConcreteACI318_25` to the factory mapping.
In `materials/reinforcement/__init__.py`, add `ReinforcementACI318_25` to the factory mapping.

## 3. Whitney Stress Block Constitutive Law

**File:** `materials/constitutive_laws/_whitneyblock.py`

```python
class WhitneyBlock(ConstitutiveLaw):
    """Equivalent rectangular stress block for section integration.

    This constitutive law represents the equivalent rectangular compressive
    stress distribution used in ACI 318 and other codes (CSA A23.3, AS 3600)
    for computing nominal flexural strength.

    It is NOT a physical stress-strain relationship. It is a code-calibrated
    design idealization that produces the same resultant force and moment as
    the actual nonlinear concrete stress distribution at nominal strength.
    The specific parameters (stress intensity, depth factor, ultimate strain)
    are code-dependent and are supplied by the material class via the
    constitutive law factory pattern (e.g., ConcreteACI318_25.__whitneyblock__()).

    For integration purposes, this is modeled as a piecewise-constant
    stress-strain function. In a linear strain profile with eps_cu at the
    extreme compression fiber:
    - Strain at depth a = beta1*c corresponds to eps_cu*(1-beta1)
    - Stress = fc for strains between eps_cu*(1-beta1) and eps_cu
    - Stress = 0 for strains between 0 and eps_cu*(1-beta1)

    This representation allows both the Marin and Fiber integrators to
    consume the Whitney block without any modification to the section
    analysis pipeline.

    Args:
        fc: Stress block intensity, typically alpha1 * f'c (MPa).
        beta1: Depth factor mapping neutral axis depth c to block depth a = beta1*c.
        eps_cu: Ultimate concrete strain (default 0.003).
    """

    __materials__ = ('concrete',)
```

Methods:
- `get_stress(eps)` — returns `-fc` in the active zone, `0` elsewhere
- `get_ultimate_strain()` — returns `(-eps_cu, 0.0)`
- `__marin__(strain)` — returns coefficients for two zones (zero-stress, constant-stress)
- `__marin_tangent__(strain)` — tangent version for Marin integration

Registered in `materials/constitutive_laws/__init__.py`:
```python
CONSTITUTIVE_LAWS = {
    ...
    'whitneyblock': WhitneyBlock,
}
```

## 4. Strength Reduction Factors

**File:** `codes/aci318_25/_strength_reduction.py`

### Fixed phi values (Table 21.2.1)

```python
def phi_shear() -> float:
    """Table 21.2.1(b). Returns 0.75."""

def phi_torsion() -> float:
    """Table 21.2.1(c). Returns 0.75."""

def phi_bearing() -> float:
    """Table 21.2.1(d). Returns 0.65."""
```

### Strain-dependent phi (Table 21.2.2)

```python
def phi_flexure(
    eps_t: float,
    fy: float,
    Es: float = 200000.0,
    transverse: Literal['spiral', 'other'] = 'other',
) -> float:
    """Strength reduction factor for moment, axial force, or combined.

    ACI 318-25, Table 21.2.2.

    Classification based on net tensile strain eps_t:
    - eps_t <= eps_ty:              compression-controlled (0.75 spiral / 0.65 other)
    - eps_ty < eps_t < eps_ty+0.003: transition (linear interpolation)
    - eps_t >= eps_ty + 0.003:      tension-controlled (0.90)

    where eps_ty = fy / Es (per 21.2.2.1).
    """
```

### Section classification helper

```python
def section_classification(
    eps_t: float, fy: float, Es: float = 200000.0,
) -> Literal['tension-controlled', 'transition', 'compression-controlled']:
    """Classify section per Table 21.2.2."""
```

## 5. Flexure Module

**File:** `codes/aci318_25/_flexure.py`

### Equilibrium helpers (shared across section types)

```python
def stress_block_depth_sr(As, fy, fc, b) -> float:
    """Stress block depth a for singly-reinforced rectangular section.
    a = As * fy / (0.85 * f'c * b)"""

def stress_block_depth_dr(As, As_prime, fy, fy_prime, fc, b) -> float:
    """Stress block depth a for doubly-reinforced rectangular section.
    a = (As*fy - As'*fy') / (0.85 * f'c * b)"""

def neutral_axis_depth(a, beta1) -> float:
    """c = a / beta1. Common to all rectangular sections."""

def eps_t_from_c(c, dt, eps_cu=0.003) -> float:
    """Net tensile strain: eps_t = eps_cu * (dt - c) / c.
    Common to all section types."""

def eps_s_prime(c, d_prime, eps_cu=0.003) -> float:
    """Compression steel strain: eps_s' = eps_cu * (c - d') / c.
    Used to verify compression steel has yielded (doubly-reinforced)."""
```

### Nominal moment strength

```python
def Mn_singly_reinforced(As, fy, fc, b, d) -> float:
    """Mn = As * fy * (d - a/2) for singly-reinforced rectangular section."""

def Mn_doubly_reinforced(As, As_prime, fy, fy_prime, fc, b, d, d_prime) -> float:
    """Mn for doubly-reinforced rectangular section.
    Mn = (As*fy - As'*fy') * (d - a/2) + As'*fy' * (d - d')
    Note: Caller must verify compression steel yields via eps_s_prime()."""
```

### Reinforcement limits (member-type-specific)

```python
def As_min_slab(fy, b, h) -> float:
    """Minimum flexural reinforcement for one-way slabs.
    7.6.1.1 -> 24.4.3.2. Same as shrinkage/temperature reinforcement."""

def As_min_beam(fc, fy, bw, d) -> float:
    """Minimum flexural reinforcement for beams.
    9.6.1.2: As_min = max(3*sqrt(f'c)/fy, 200/fy) * bw * d"""

def As_max_check(eps_t, fy, Es=200000.0) -> bool:
    """Check tension-controlled: eps_t >= eps_ty + 0.003.
    Required for slabs (7.3.3.1) and beams (9.3.3.1)."""
```

### Design helpers

```python
def As_required(Mu, phi, fy, fc, b, d) -> float:
    """Required As for singly-reinforced section given factored Mu.
    Quadratic solution: Mu = phi * As * fy * (d - As*fy / (1.7*f'c*b))"""
```

## 6. Shear Module

**File:** `codes/aci318_25/_shear.py`

### Size effect

```python
def lambda_s(d: float) -> float:
    """Size effect modification factor. Eq. 22.5.5.1.3.
    lambda_s = 2 / (1 + d/10) <= 1.0
    Note: d is in inches in the code. This function accepts mm and converts internally."""
```

### Concrete shear contribution

```python
def Vc_detailed(
    fc, bw, d, rho_w, Nu=0.0, Ag=0.0, lambda_concrete=1.0,
    Av_provided=0.0, Av_min=0.0,
) -> float:
    """Concrete shear strength for nonprestressed members. Table 22.5.5.1.

    If Av >= Av_min:
        Vc = [8*lambda*(rho_w)^(1/3)*sqrt(f'c) + Nu/(6*Ag)] * bw * d   (b)
    If Av < Av_min:
        Vc = [8*lambda_s*lambda*(rho_w)^(1/3)*sqrt(f'c) + Nu/(6*Ag)] * bw * d   (c)

    Limits per 22.5.5.1.1:
        Vc <= 5*lambda*sqrt(f'c)*bw*d
        Vc >= lambda*sqrt(f'c)*bw*d  (unless net axial tension)
    Per 22.5.5.1.2: Nu/(6*Ag) <= 0.05*f'c
    Per 22.5.3.1: sqrt(f'c) <= 8.3 MPa (100 psi equivalent)
    Per 22.5.3.3: fy, fyt <= 420 MPa for shear calcs
    """

def Vc_simplified(fc, bw, d, Nu=0.0, Ag=0.0, lambda_concrete=1.0) -> float:
    """Simplified Vc. Table 22.5.5.1(a).
    Vc = [2*lambda*sqrt(f'c) + Nu/(6*Ag)] * bw * d
    Only valid when Av >= Av_min."""
```

### Steel shear contribution

```python
def Vs(Av, fyt, d, s) -> float:
    """Eq. 22.5.8.5.3. Vs = Av * fyt * d / s"""

def Vn(Vc, Vs) -> float:
    """Nominal shear strength. Vn = Vc + Vs."""
```

### Checks and limits

```python
def check_cross_section(Vu, phi, Vc, fc, bw, d) -> bool:
    """Eq. 22.5.1.2: Vu <= phi * (Vc + 8*sqrt(f'c)*bw*d)"""

def Av_min_per_s(fc, bw, fyt) -> float:
    """Minimum shear reinforcement. 9.6.3.4 / 7.6.3.3.
    Av_min/s = max(0.062*sqrt(f'c), 0.35) * bw / fyt"""

def shear_reinforcement_required(Vu, phi_Vc) -> bool:
    """For slabs (7.6.3.1): required when Vu > phi*Vc."""

def max_stirrup_spacing(d, Vs, fc, bw) -> float:
    """9.7.6.2.2. d/2 or d/4 depending on Vs level."""
```

## 7. One-Way Slab Module

**File:** `codes/aci318_25/_one_way_slab.py`

Slab-specific rules from Ch. 7 that reference the shared flexure/shear functions.

```python
def min_thickness(span, support_condition, fy=420.0, lightweight=False, wc=2320.0) -> float:
    """Table 7.3.1.1. L/20, L/24, L/28, L/10 with adjustments for fy and lightweight."""

def As_shrinkage_temperature(fy, b, h) -> float:
    """24.4.3.2. Gr 60: 0.0018*b*h. Also the minimum flexural reinforcement (7.6.1.1)."""

def max_bar_spacing_flexure(h) -> float:
    """7.7.2.3: min(3*h, 450 mm)."""

def max_bar_spacing_shrinkage(h) -> float:
    """7.7.6.2.1: min(5*h, 450 mm)."""

def shear_critical_section_offset(d) -> float:
    """7.4.3.2: d from face of support for nonprestressed slabs."""
```

### Future member modules (same pattern)

```
codes/aci318_25/
├── _beam.py             # Ch. 9 (future)
├── _column.py           # Ch. 10 (future)
├── _wall.py             # Ch. 11 (future)
├── _two_way_slab.py     # Ch. 8 (future)
└── _foundation.py       # Ch. 13 (future)
```

Each member module applies the correct limits, minimums, and detailing from its ACI
chapter, calling the shared flexure/shear/phi functions underneath.

## 8. Files Modified vs. Created

### Existing files modified (minimal, registration only)

| File | Change |
|------|--------|
| `codes/__init__.py` | Add `aci318_25` to imports and `_DESIGN_CODES` dict |
| `materials/concrete/__init__.py` | Add `ConcreteACI318_25` to factory mapping |
| `materials/reinforcement/__init__.py` | Add `ReinforcementACI318_25` to factory mapping |
| `materials/constitutive_laws/__init__.py` | Add `WhitneyBlock` to `CONSTITUTIVE_LAWS` dict |

### New files created

| File | Purpose |
|------|---------|
| `codes/aci318_25/__init__.py` | Module metadata and public API |
| `codes/aci318_25/_concrete_material_properties.py` | Ch. 19 material property functions |
| `codes/aci318_25/_reinforcement_material_properties.py` | Ch. 20 material property functions |
| `codes/aci318_25/_strength_reduction.py` | Ch. 21 phi factors |
| `codes/aci318_25/_flexure.py` | Ch. 22.2-22.3 flexural strength |
| `codes/aci318_25/_shear.py` | Ch. 22.5 one-way shear strength |
| `codes/aci318_25/_one_way_slab.py` | Ch. 7 slab-specific rules |
| `codes/aci318_25/_units.py` | Unit conversion constants (PSI_TO_MPA, IN_TO_MM, FT_TO_MM, etc.) |
| `materials/concrete/_concreteACI318_25.py` | Concrete material class |
| `materials/reinforcement/_reinforcementACI318_25.py` | Reinforcement material class |
| `materials/constitutive_laws/_whitneyblock.py` | Whitney stress block constitutive law |

### Untouched

- `geometry/` — all geometry modules
- `sections/` — `BeamSection`, `BeamSectionCalculator`, integrators
- `core/` — base classes
- All existing code modules (`ec2_2004/`, `ec2_2023/`, `mc2010/`, `mc2020/`)
- All existing constitutive laws
- All existing material base classes

## 9. Example Problem — One-Way Slab Validation

**Parameters:**
- f'c = 4000 psi = 27.58 MPa
- fy = 60 ksi = 420 MPa
- Span: 20 ft = 6096 mm (clear span), one end continuous
- Loads: self-weight + 20 psf SDL + 100 psf LL (UDL), plus 5 kip point load at midspan
- Width: design per 1 ft strip (b = 305 mm)

### Path A — Closed-form

```python
from structuralcodes.codes import aci318_25

# Thickness
h = aci318_25.min_thickness(span=6096, support_condition='one_end_continuous')  # L/24

# Effective depth (#5 bars, 3/4" cover)
d = h - 19 - 16/2  # ~227 mm

# Material properties
beta1 = aci318_25.beta1(27.58)  # 0.85

# Flexure (after computing Mu from load combinations)
As_req = aci318_25.As_required(Mu=Mu, phi=0.9, fy=420, fc=27.58, b=305, d=227)
As_min = aci318_25.As_min_slab(fy=420, b=305, h=254)

# Tension-controlled check
a = aci318_25.stress_block_depth_sr(As=As_req, fy=420, fc=27.58, b=305)
c = aci318_25.neutral_axis_depth(a=a, beta1=0.85)
eps_t = aci318_25.eps_t_from_c(c=c, dt=227)
phi = aci318_25.phi_flexure(eps_t=eps_t, fy=420)

# Shear
Vc = aci318_25.Vc_detailed(fc=27.58, bw=305, d=227, rho_w=As_req/(305*227))
assert Vu <= aci318_25.phi_shear() * Vc  # No stirrups needed
```

### Path B — Section integrator

```python
from structuralcodes.materials.concrete import ConcreteACI318_25
from structuralcodes.materials.reinforcement import ReinforcementACI318_25
from structuralcodes.geometry import SurfaceGeometry, PointGeometry, CompoundGeometry
from structuralcodes.sections import BeamSection
from structuralcodes.codes import aci318_25

concrete = ConcreteACI318_25(fck=27.58, constitutive_law='parabolarectangle')
steel = ReinforcementACI318_25(fyk=420, Es=200000, ftk=550, epsuk=0.05,
                                constitutive_law='elasticperfectlyplastic')

poly = Polygon([(0, 0), (305, 0), (305, 254), (0, 254)])
surf = SurfaceGeometry(poly, concrete)
bar1 = PointGeometry(point=(100, 27), diameter=16, material=steel)
bar2 = PointGeometry(point=(205, 27), diameter=16, material=steel)
section_geo = CompoundGeometry([surf], [bar1, bar2])

section = BeamSection(section_geo, integrator='marin')
calc = section.section_calculator
strain_profile = calc.find_equilibrium_fixed_pivot(geom=section.geometry, n=0, yielding=True)
N, My, Mz, data = calc.integrator.integrate_strain_response_on_geometry(
    geometry=section.geometry, strain=strain_profile
)

Mn_integrator = abs(My)
# Extract eps_t from the strain profile at the extreme tension steel location
eps_0, kappa = strain_profile[0], strain_profile[1]
eps_t = eps_0 + kappa * (254 - 27)  # strain at bottom reinforcement level
phi = aci318_25.phi_flexure(eps_t=eps_t, fy=420)
phi_Mn = phi * Mn_integrator
# Cross-check: should agree with Path A within ~2-5%
```

## 10. Test Structure

```
tests/test_aci318_25/
├── __init__.py
├── test_concrete_material_properties.py      # Ec, fr, beta1, lambda_factor, eps_cu, alpha1, fct
├── test_reinforcement_material_properties.py  # Es, fy_design, epsyd, grade lookup
├── test_concrete_aci318_25.py                # Material class, constitutive law creation, factory
├── test_reinforcement_aci318_25.py           # Material class, from_grade, factory
├── test_strength_reduction.py                # phi_flexure, phi_shear, transition zone, classification
├── test_flexure.py                           # Mn (SR/DR), As_required, As_min, eps_t, eps_s_prime
├── test_shear.py                             # Vc_detailed, Vc_simplified, Vs, Vn, lambda_s, limits
├── test_one_way_slab.py                      # min_thickness, shrinkage steel, spacing limits
├── test_whitneyblock.py                      # get_stress, get_ultimate_strain, __marin__, integration
└── test_one_way_slab_example.py              # End-to-end: both paths, cross-check results
```

## References

- ACI 318-25: Building Code Requirements and Commentary for Structural Concrete
  (available at `\\KPL\VA-Prj$\1\03\00764\01\A\Data\5_References\Technical Resources\ACI\`)
- structuralcodes repository: https://github.com/fib-international/structuralcodes
- Issue #187: https://github.com/fib-international/structuralcodes/issues/187
- PR #343 (gabe-kafka, ACI 318-19 materials): https://github.com/fib-international/structuralcodes/pull/343
