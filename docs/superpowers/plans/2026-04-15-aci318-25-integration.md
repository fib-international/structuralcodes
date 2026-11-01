# ACI 318-25 Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add ACI 318-25 one-way slab flexural and shear design to structuralcodes, with both closed-form equations and section-integrator paths.

**Architecture:** Thin material adapter classes (`ConcreteACI318_25`, `ReinforcementACI318_25`) plug into the existing geometry/section/integrator pipeline unchanged. ACI design intelligence (flexure, shear, phi factors) lives in pure functions under `codes/aci318_25/`. A `WhitneyBlock` constitutive law is added to `materials/constitutive_laws/` for integrator-driven analysis.

**Tech Stack:** Python 3.10+, NumPy, Shapely, pytest. SI units throughout (MPa, mm, N).

**Spec:** `docs/superpowers/specs/2026-04-15-aci318-25-integration-design.md`

---

## File Map

### New files

| File | Responsibility |
|------|---------------|
| `structuralcodes/codes/aci318_25/__init__.py` | Module metadata, public API exports |
| `structuralcodes/codes/aci318_25/_concrete_material_properties.py` | Ch. 19: `Ec`, `fr`, `beta1`, `eps_cu`, `alpha1`, `fct`, `lambda_factor` |
| `structuralcodes/codes/aci318_25/_reinforcement_material_properties.py` | Ch. 20: `Es`, `fy_design`, `epsyd`, `reinforcement_grade_props` |
| `structuralcodes/codes/aci318_25/_strength_reduction.py` | Ch. 21: `phi_flexure`, `phi_shear`, `section_classification` |
| `structuralcodes/codes/aci318_25/_flexure.py` | Ch. 22.2-22.3: `Mn_singly_reinforced`, `Mn_doubly_reinforced`, `As_required`, helpers |
| `structuralcodes/codes/aci318_25/_shear.py` | Ch. 22.5: `Vc_detailed`, `Vc_simplified`, `Vs`, `Vn`, limits |
| `structuralcodes/codes/aci318_25/_one_way_slab.py` | Ch. 7: `min_thickness`, `As_shrinkage_temperature`, spacing limits |
| `structuralcodes/codes/aci318_25/_units.py` | Conversion constants: `PSI_TO_MPA`, `IN_TO_MM`, etc. |
| `structuralcodes/materials/concrete/_concreteACI318_25.py` | `ConcreteACI318_25` material class |
| `structuralcodes/materials/reinforcement/_reinforcementACI318_25.py` | `ReinforcementACI318_25` material class |
| `structuralcodes/materials/constitutive_laws/_whitneyblock.py` | `WhitneyBlock` constitutive law |
| `tests/test_aci318_25/__init__.py` | Test package |
| `tests/test_aci318_25/test_concrete_material_properties.py` | Tests for Ch. 19 functions |
| `tests/test_aci318_25/test_reinforcement_material_properties.py` | Tests for Ch. 20 functions |
| `tests/test_aci318_25/test_strength_reduction.py` | Tests for phi factors |
| `tests/test_aci318_25/test_flexure.py` | Tests for flexure functions |
| `tests/test_aci318_25/test_shear.py` | Tests for shear functions |
| `tests/test_aci318_25/test_one_way_slab.py` | Tests for slab rules |
| `tests/test_aci318_25/test_concrete_aci318_25.py` | Tests for concrete material class |
| `tests/test_aci318_25/test_reinforcement_aci318_25.py` | Tests for reinforcement material class |
| `tests/test_aci318_25/test_whitneyblock.py` | Tests for Whitney block constitutive law |
| `tests/test_aci318_25/test_one_way_slab_example.py` | End-to-end validation (both paths) |

### Modified files (registration only)

| File | Change |
|------|--------|
| `structuralcodes/codes/__init__.py` | Add `aci318_25` import and registry entry |
| `structuralcodes/materials/concrete/__init__.py` | Add `ConcreteACI318_25` import and `CONCRETES` entry |
| `structuralcodes/materials/reinforcement/__init__.py` | Add `ReinforcementACI318_25` import and `REINFORCEMENTS` entry |
| `structuralcodes/materials/constitutive_laws/__init__.py` | Add `WhitneyBlock` import and `CONSTITUTIVE_LAWS` entry |

---

## Task 1: Unit Conversion Constants and Code Module Skeleton

**Files:**
- Create: `structuralcodes/codes/aci318_25/__init__.py`
- Create: `structuralcodes/codes/aci318_25/_units.py`
- Modify: `structuralcodes/codes/__init__.py`

- [ ] **Step 1: Create the `_units.py` module**

```python
# structuralcodes/codes/aci318_25/_units.py
"""Unit conversion constants for ACI 318-25.

The structuralcodes library uses SI units internally (MPa, mm, N, kg/m3).
ACI 318 is published in US customary units. These constants allow users
to convert between systems at the API boundary.
"""

# Stress
PSI_TO_MPA = 0.00689476
KSI_TO_MPA = 6.89476
MPA_TO_PSI = 145.038
MPA_TO_KSI = 0.145038

# Length
IN_TO_MM = 25.4
FT_TO_MM = 304.8
MM_TO_IN = 1.0 / 25.4
MM_TO_FT = 1.0 / 304.8

# Force
LBF_TO_N = 4.44822
KIP_TO_N = 4448.22
N_TO_LBF = 1.0 / 4.44822
N_TO_KIP = 1.0 / 4448.22

# Distributed load
PSF_TO_PA = 47.8803
PSF_TO_KPA = 0.0478803

# Density
PCF_TO_KGM3 = 16.0185
KGM3_TO_PCF = 1.0 / 16.0185
```

- [ ] **Step 2: Create the code module `__init__.py`**

```python
# structuralcodes/codes/aci318_25/__init__.py
"""ACI 318-25: Building Code Requirements for Structural Concrete."""

import typing as t

from ._units import (
    FT_TO_MM,
    IN_TO_MM,
    KIP_TO_N,
    KSI_TO_MPA,
    PSI_TO_MPA,
)

__all__: t.List[str] = [
    'PSI_TO_MPA',
    'KSI_TO_MPA',
    'IN_TO_MM',
    'FT_TO_MM',
    'KIP_TO_N',
]

__title__: str = 'ACI 318-25'
__year__: str = '2025'
__materials__: t.Tuple[str, ...] = ('concrete', 'reinforcement')
```

- [ ] **Step 3: Register in the design code registry**

In `structuralcodes/codes/__init__.py`, add the import and registry entry.

Change the import line from:
```python
from . import ec2_2004, ec2_2023, mc2010, mc2020
```
to:
```python
from . import aci318_25, ec2_2004, ec2_2023, mc2010, mc2020
```

Add `'aci318_25'` to `__all__`:
```python
__all__ = [
    'aci318_25',
    'mc2010',
    'mc2020',
    'ec2_2023',
    'ec2_2004',
    'set_design_code',
    'get_design_codes',
    'set_national_annex',
]
```

Add to `_DESIGN_CODES`:
```python
_DESIGN_CODES = {
    'aci318_25': aci318_25,
    'mc2010': mc2010,
    'mc2020': mc2020,
    'ec2_2004': ec2_2004,
    'ec2_2023': ec2_2023,
}
```

- [ ] **Step 4: Verify the module loads**

Run: `python -c "import structuralcodes; print(structuralcodes.get_design_codes())"`

Expected output should include `'aci318_25'` in the list.

- [ ] **Step 5: Commit**

```bash
git add structuralcodes/codes/aci318_25/__init__.py structuralcodes/codes/aci318_25/_units.py structuralcodes/codes/__init__.py
git commit -m "feat(aci318_25): add code module skeleton and unit conversion constants"
```

---

## Task 2: Concrete Material Property Functions

**Files:**
- Create: `structuralcodes/codes/aci318_25/_concrete_material_properties.py`
- Create: `tests/test_aci318_25/__init__.py`
- Create: `tests/test_aci318_25/test_concrete_material_properties.py`
- Modify: `structuralcodes/codes/aci318_25/__init__.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_aci318_25/__init__.py
"""Collection of tests for ACI 318-25."""
```

```python
# tests/test_aci318_25/test_concrete_material_properties.py
"""Tests for concrete material properties of ACI 318-25."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _concrete_material_properties as cmp


class TestEc:
    """Tests for modulus of elasticity (Table 19.2.2.1)."""

    def test_normalweight_4000psi(self):
        """Ec for f'c = 27.58 MPa (4000 psi), wc = 2320 kg/m3."""
        expected = 2320**1.5 * 0.043 * math.sqrt(27.58)
        assert math.isclose(cmp.Ec(27.58), expected, rel_tol=1e-6)

    def test_normalweight_28mpa(self):
        """Ec for f'c = 28 MPa."""
        expected = 2320**1.5 * 0.043 * math.sqrt(28)
        assert math.isclose(cmp.Ec(28), expected, rel_tol=1e-6)

    def test_custom_unit_weight(self):
        """Ec with lightweight concrete wc = 1800 kg/m3."""
        expected = 1800**1.5 * 0.043 * math.sqrt(28)
        assert math.isclose(cmp.Ec(28, wc=1800), expected, rel_tol=1e-6)

    def test_invalid_fc_raises(self):
        with pytest.raises(ValueError):
            cmp.Ec(-1)

    def test_invalid_wc_raises(self):
        with pytest.raises(ValueError):
            cmp.Ec(28, wc=1000)


class TestFr:
    """Tests for modulus of rupture (Eq. 19.2.3.1)."""

    def test_normalweight(self):
        expected = 0.62 * math.sqrt(27.58)
        assert math.isclose(cmp.fr(27.58), expected, rel_tol=1e-6)

    def test_lightweight(self):
        expected = 0.62 * 0.75 * math.sqrt(27.58)
        assert math.isclose(cmp.fr(27.58, lambda_s=0.75), expected, rel_tol=1e-6)

    def test_invalid_fc_raises(self):
        with pytest.raises(ValueError):
            cmp.fr(-1)

    def test_invalid_lambda_raises(self):
        with pytest.raises(ValueError):
            cmp.fr(28, lambda_s=1.5)


class TestBeta1:
    """Tests for Whitney stress block depth factor (Table 22.2.2.4.3)."""

    @pytest.mark.parametrize('fc, expected', [
        (21, 0.85),       # Below 28 MPa (4000 psi)
        (27.58, 0.85),    # At 28 MPa (4000 psi)
        (34.47, 0.8036),  # 5000 psi = 34.47 MPa: 0.85 - 0.05*(34.47-28)/7
        (55.16, 0.65),    # 8000 psi = 55.16 MPa
        (68.95, 0.65),    # Above 55 MPa
    ])
    def test_beta1_values(self, fc, expected):
        assert math.isclose(cmp.beta1(fc), expected, rel_tol=1e-3)

    def test_invalid_fc_raises(self):
        with pytest.raises(ValueError):
            cmp.beta1(-1)


class TestEpsCu:
    """Tests for ultimate concrete strain (Sec. 22.2.2.1)."""

    def test_value(self):
        assert cmp.eps_cu() == 0.003


class TestAlpha1:
    """Tests for stress block intensity (Sec. 22.2.2.4.1)."""

    def test_value(self):
        assert cmp.alpha1() == 0.85


class TestFct:
    """Tests for splitting tensile strength (Sec. 19.2.4.3)."""

    def test_normalweight(self):
        expected = 0.56 * math.sqrt(27.58)
        assert math.isclose(cmp.fct(27.58), expected, rel_tol=1e-6)

    def test_invalid_fc_raises(self):
        with pytest.raises(ValueError):
            cmp.fct(-1)


class TestLambdaFactor:
    """Tests for lightweight modification factor (Table 19.2.4.2)."""

    @pytest.mark.parametrize('concrete_type, expected', [
        ('normalweight', 1.0),
        ('sand-lightweight', 0.85),
        ('all-lightweight', 0.75),
    ])
    def test_known_types(self, concrete_type, expected):
        assert cmp.lambda_factor(concrete_type) == expected

    def test_invalid_type_raises(self):
        with pytest.raises(ValueError):
            cmp.lambda_factor('unknown')
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_aci318_25/test_concrete_material_properties.py -v`

Expected: ERRORS — module `_concrete_material_properties` does not exist.

- [ ] **Step 3: Implement concrete material property functions**

```python
# structuralcodes/codes/aci318_25/_concrete_material_properties.py
"""Concrete material properties according to ACI 318-25, Chapter 19."""

from __future__ import annotations

import math
import typing as t

LAMBDA_FACTORS = {
    'normalweight': 1.0,
    'sand-lightweight': 0.85,
    'all-lightweight': 0.75,
}


def Ec(fc: float, wc: float = 2320.0) -> float:
    """Modulus of elasticity of concrete.

    ACI 318-25, Table 19.2.2.1.

    Args:
        fc: Specified compressive strength f'c in MPa.
        wc: Unit weight of concrete in kg/m3 (default 2320, normalweight).

    Returns:
        Modulus of elasticity in MPa.

    Raises:
        ValueError: If fc is not positive.
        ValueError: If wc is outside 1440-2560 kg/m3.
    """
    if fc <= 0:
        raise ValueError(f'fc={fc} must be positive')
    if wc < 1440 or wc > 2560:
        raise ValueError(f'wc={wc} must be between 1440 and 2560 kg/m3')
    return wc**1.5 * 0.043 * math.sqrt(fc)


def fr(fc: float, lambda_s: float = 1.0) -> float:
    """Modulus of rupture of concrete.

    ACI 318-25, Eq. 19.2.3.1.

    Args:
        fc: Specified compressive strength f'c in MPa.
        lambda_s: Lightweight concrete modification factor (default 1.0).

    Returns:
        Modulus of rupture in MPa.

    Raises:
        ValueError: If fc is not positive.
        ValueError: If lambda_s is not in (0, 1].
    """
    if fc <= 0:
        raise ValueError(f'fc={fc} must be positive')
    if lambda_s <= 0 or lambda_s > 1.0:
        raise ValueError(f'lambda_s={lambda_s} must be in the range (0, 1]')
    return 0.62 * lambda_s * math.sqrt(fc)


def beta1(fc: float) -> float:
    """Whitney stress block depth factor.

    ACI 318-25, Table 22.2.2.4.3.

    Args:
        fc: Specified compressive strength f'c in MPa.

    Returns:
        Stress block depth factor (dimensionless).

    Raises:
        ValueError: If fc is not positive.
    """
    if fc <= 0:
        raise ValueError(f'fc={fc} must be positive')
    if fc <= 28:
        return 0.85
    if fc >= 55:
        return 0.65
    return 0.85 - 0.05 * (fc - 28) / 7


def eps_cu() -> float:
    """Maximum usable strain at extreme concrete compression fiber.

    ACI 318-25, Sec. 22.2.2.1.

    Returns:
        Ultimate concrete strain (dimensionless).
    """
    return 0.003


def alpha1() -> float:
    """Ratio of equivalent rectangular stress block intensity.

    ACI 318-25, Sec. 22.2.2.4.1.

    Returns:
        Stress block intensity factor (dimensionless).
    """
    return 0.85


def fct(fc: float, lambda_s: float = 1.0) -> float:
    """Approximate splitting tensile strength of concrete.

    ACI 318-25, Sec. 19.2.4.3.

    Args:
        fc: Specified compressive strength f'c in MPa.
        lambda_s: Lightweight concrete modification factor (default 1.0).

    Returns:
        Splitting tensile strength in MPa.

    Raises:
        ValueError: If fc is not positive.
    """
    if fc <= 0:
        raise ValueError(f'fc={fc} must be positive')
    return 0.56 * lambda_s * math.sqrt(fc)


def lambda_factor(
    concrete_type: t.Literal[
        'normalweight', 'sand-lightweight', 'all-lightweight'
    ],
) -> float:
    """Lightweight concrete modification factor.

    ACI 318-25, Table 19.2.4.2.

    Args:
        concrete_type: One of 'normalweight', 'sand-lightweight',
            or 'all-lightweight'.

    Returns:
        Lightweight modification factor (dimensionless).

    Raises:
        ValueError: If concrete_type is not recognized.
    """
    result = LAMBDA_FACTORS.get(concrete_type.lower())
    if result is None:
        raise ValueError(
            f'Unknown concrete type: {concrete_type}. '
            f'Valid types: {list(LAMBDA_FACTORS.keys())}'
        )
    return result
```

- [ ] **Step 4: Update `__init__.py` to export these functions**

Replace the contents of `structuralcodes/codes/aci318_25/__init__.py` with:

```python
# structuralcodes/codes/aci318_25/__init__.py
"""ACI 318-25: Building Code Requirements for Structural Concrete."""

import typing as t

from ._concrete_material_properties import (
    Ec,
    alpha1,
    beta1,
    eps_cu,
    fct,
    fr,
    lambda_factor,
)
from ._units import (
    FT_TO_MM,
    IN_TO_MM,
    KIP_TO_N,
    KSI_TO_MPA,
    PSI_TO_MPA,
)

__all__: t.List[str] = [
    'Ec',
    'alpha1',
    'beta1',
    'eps_cu',
    'fct',
    'fr',
    'lambda_factor',
    'PSI_TO_MPA',
    'KSI_TO_MPA',
    'IN_TO_MM',
    'FT_TO_MM',
    'KIP_TO_N',
]

__title__: str = 'ACI 318-25'
__year__: str = '2025'
__materials__: t.Tuple[str, ...] = ('concrete', 'reinforcement')
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_aci318_25/test_concrete_material_properties.py -v`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add structuralcodes/codes/aci318_25/_concrete_material_properties.py structuralcodes/codes/aci318_25/__init__.py tests/test_aci318_25/
git commit -m "feat(aci318_25): add concrete material property functions (Ch. 19)"
```

---

## Task 3: Reinforcement Material Property Functions

**Files:**
- Create: `structuralcodes/codes/aci318_25/_reinforcement_material_properties.py`
- Create: `tests/test_aci318_25/test_reinforcement_material_properties.py`
- Modify: `structuralcodes/codes/aci318_25/__init__.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_aci318_25/test_reinforcement_material_properties.py
"""Tests for reinforcement material properties of ACI 318-25."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _reinforcement_material_properties as rmp


class TestEs:
    def test_value(self):
        assert rmp.Es() == 200000.0


class TestFyDesign:
    def test_default_phi(self):
        assert math.isclose(rmp.fy_design(420), 420)

    def test_with_phi(self):
        assert math.isclose(rmp.fy_design(420, phi=0.9), 378)

    def test_invalid_fy_raises(self):
        with pytest.raises(ValueError):
            rmp.fy_design(-1)

    def test_invalid_phi_raises(self):
        with pytest.raises(ValueError):
            rmp.fy_design(420, phi=1.5)


class TestEpsyd:
    @pytest.mark.parametrize('fy, expected', [
        (420, 420 / 200000),
        (280, 280 / 200000),
        (550, 550 / 200000),
    ])
    def test_yield_strain(self, fy, expected):
        assert math.isclose(rmp.epsyd(fy), expected, rel_tol=1e-6)

    def test_invalid_fy_raises(self):
        with pytest.raises(ValueError):
            rmp.epsyd(-1)


class TestReinforcementGradeProps:
    @pytest.mark.parametrize('grade, exp_fy, exp_fu', [
        ('40', 280.0, 420.0),
        ('60', 420.0, 550.0),
        ('80', 550.0, 690.0),
        ('100', 690.0, 860.0),
    ])
    def test_known_grades(self, grade, exp_fy, exp_fu):
        props = rmp.reinforcement_grade_props(grade)
        assert math.isclose(props['fy'], exp_fy)
        assert math.isclose(props['fu'], exp_fu)

    def test_invalid_grade_raises(self):
        with pytest.raises(ValueError):
            rmp.reinforcement_grade_props('999')
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_aci318_25/test_reinforcement_material_properties.py -v`

Expected: ERRORS — module does not exist.

- [ ] **Step 3: Implement reinforcement material property functions**

```python
# structuralcodes/codes/aci318_25/_reinforcement_material_properties.py
"""Reinforcement material properties according to ACI 318-25, Chapter 20."""

from __future__ import annotations

import typing as t

REINFORCEMENT_GRADES: t.Dict[str, t.Dict[str, float]] = {
    '40': {'fy': 280.0, 'fu': 420.0},
    '60': {'fy': 420.0, 'fu': 550.0},
    '80': {'fy': 550.0, 'fu': 690.0},
    '100': {'fy': 690.0, 'fu': 860.0},
}


def Es() -> float:
    """Modulus of elasticity of reinforcement.

    ACI 318-25, Sec. 20.2.2.2.

    Returns:
        Modulus of elasticity in MPa.
    """
    return 200000.0


def fy_design(fy: float, phi: float = 1.0) -> float:
    """Design yield strength of reinforcement.

    ACI 318-25 applies strength reduction factors (phi) at the member
    capacity level, not the material level. Default phi=1.0 returns
    the unreduced yield strength.

    Args:
        fy: Specified yield strength in MPa.
        phi: Optional strength reduction factor (default 1.0).

    Returns:
        Design yield strength in MPa.

    Raises:
        ValueError: If fy is not positive.
        ValueError: If phi is not in (0, 1].
    """
    if fy <= 0:
        raise ValueError(f'fy={fy} must be positive')
    if phi <= 0 or phi > 1.0:
        raise ValueError(f'phi={phi} must be in the range (0, 1]')
    return phi * fy


def epsyd(fy: float, _Es: float = 200000.0) -> float:
    """Yield strain of reinforcement.

    Args:
        fy: Specified yield strength in MPa.
        _Es: Modulus of elasticity in MPa (default 200000).

    Returns:
        Yield strain (dimensionless).

    Raises:
        ValueError: If fy is not positive.
    """
    if fy <= 0:
        raise ValueError(f'fy={fy} must be positive')
    return fy / _Es


def reinforcement_grade_props(grade: str) -> t.Dict[str, float]:
    """Look up ASTM A615 reinforcement grade properties.

    ACI 318-25, Table 20.2.2.4a.

    Args:
        grade: ASTM grade as string ('40', '60', '80', '100').

    Returns:
        Dict with 'fy' (MPa) and 'fu' (MPa).

    Raises:
        ValueError: If grade is not recognized.
    """
    props = REINFORCEMENT_GRADES.get(grade)
    if props is None:
        raise ValueError(
            f'Unknown grade: {grade}. '
            f'Valid grades: {list(REINFORCEMENT_GRADES.keys())}'
        )
    return dict(props)
```

- [ ] **Step 4: Update `__init__.py` to export these functions**

Add to `structuralcodes/codes/aci318_25/__init__.py` the new imports:

```python
from ._reinforcement_material_properties import (
    Es,
    epsyd,
    fy_design,
    reinforcement_grade_props,
)
```

And add to `__all__`:
```python
'Es',
'epsyd',
'fy_design',
'reinforcement_grade_props',
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_aci318_25/test_reinforcement_material_properties.py -v`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add structuralcodes/codes/aci318_25/_reinforcement_material_properties.py structuralcodes/codes/aci318_25/__init__.py tests/test_aci318_25/test_reinforcement_material_properties.py
git commit -m "feat(aci318_25): add reinforcement material property functions (Ch. 20)"
```

---

## Task 4: Strength Reduction Factors

**Files:**
- Create: `structuralcodes/codes/aci318_25/_strength_reduction.py`
- Create: `tests/test_aci318_25/test_strength_reduction.py`
- Modify: `structuralcodes/codes/aci318_25/__init__.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_aci318_25/test_strength_reduction.py
"""Tests for strength reduction factors of ACI 318-25, Chapter 21."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _strength_reduction as sr


class TestPhiShear:
    def test_value(self):
        assert sr.phi_shear() == 0.75


class TestPhiTorsion:
    def test_value(self):
        assert sr.phi_torsion() == 0.75


class TestPhiBearing:
    def test_value(self):
        assert sr.phi_bearing() == 0.65


class TestPhiFlexure:
    """Tests for Table 21.2.2."""

    def test_tension_controlled_gr60(self):
        """eps_t = 0.005 >> eps_ty + 0.003 = 0.0051 for Gr 60."""
        # eps_ty = 420/200000 = 0.0021, limit = 0.0051
        # eps_t = 0.010 is well above limit
        assert sr.phi_flexure(eps_t=0.010, fy=420) == 0.90

    def test_tension_controlled_at_limit_gr60(self):
        """eps_t exactly at eps_ty + 0.003."""
        eps_ty = 420 / 200000  # 0.0021
        assert sr.phi_flexure(eps_t=eps_ty + 0.003, fy=420) == 0.90

    def test_compression_controlled_gr60(self):
        """eps_t <= eps_ty, other transverse."""
        eps_ty = 420 / 200000
        assert sr.phi_flexure(eps_t=eps_ty, fy=420, transverse='other') == 0.65

    def test_compression_controlled_spiral(self):
        """eps_t <= eps_ty, spiral transverse."""
        eps_ty = 420 / 200000
        assert sr.phi_flexure(eps_t=eps_ty, fy=420, transverse='spiral') == 0.75

    def test_transition_zone_midpoint_other(self):
        """Midpoint of transition zone, other transverse.
        phi = 0.65 + 0.25 * (eps_t - eps_ty) / 0.003
        At midpoint: eps_t = eps_ty + 0.0015
        phi = 0.65 + 0.25 * 0.0015 / 0.003 = 0.65 + 0.125 = 0.775
        """
        eps_ty = 420 / 200000
        phi = sr.phi_flexure(eps_t=eps_ty + 0.0015, fy=420, transverse='other')
        assert math.isclose(phi, 0.775, rel_tol=1e-6)

    def test_transition_zone_midpoint_spiral(self):
        """phi = 0.75 + 0.15 * 0.0015 / 0.003 = 0.75 + 0.075 = 0.825"""
        eps_ty = 420 / 200000
        phi = sr.phi_flexure(eps_t=eps_ty + 0.0015, fy=420, transverse='spiral')
        assert math.isclose(phi, 0.825, rel_tol=1e-6)

    def test_gr80_tension_controlled(self):
        """Gr 80: eps_ty = 550/200000 = 0.00275, limit = 0.00575."""
        assert sr.phi_flexure(eps_t=0.010, fy=550) == 0.90

    def test_gr80_compression_controlled(self):
        eps_ty = 550 / 200000
        assert sr.phi_flexure(eps_t=eps_ty, fy=550, transverse='other') == 0.65


class TestSectionClassification:
    def test_tension_controlled(self):
        eps_ty = 420 / 200000
        assert sr.section_classification(eps_ty + 0.003, fy=420) == 'tension-controlled'

    def test_compression_controlled(self):
        eps_ty = 420 / 200000
        assert sr.section_classification(eps_ty, fy=420) == 'compression-controlled'

    def test_transition(self):
        eps_ty = 420 / 200000
        assert sr.section_classification(eps_ty + 0.001, fy=420) == 'transition'
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_aci318_25/test_strength_reduction.py -v`

Expected: ERRORS — module does not exist.

- [ ] **Step 3: Implement strength reduction functions**

```python
# structuralcodes/codes/aci318_25/_strength_reduction.py
"""Strength reduction factors according to ACI 318-25, Chapter 21."""

from __future__ import annotations

import typing as t


def phi_shear() -> float:
    """Strength reduction factor for shear.

    ACI 318-25, Table 21.2.1(b).

    Returns:
        0.75
    """
    return 0.75


def phi_torsion() -> float:
    """Strength reduction factor for torsion.

    ACI 318-25, Table 21.2.1(c).

    Returns:
        0.75
    """
    return 0.75


def phi_bearing() -> float:
    """Strength reduction factor for bearing.

    ACI 318-25, Table 21.2.1(d).

    Returns:
        0.65
    """
    return 0.65


def phi_flexure(
    eps_t: float,
    fy: float,
    Es: float = 200000.0,
    transverse: t.Literal['spiral', 'other'] = 'other',
) -> float:
    """Strength reduction factor for moment, axial force, or combined.

    ACI 318-25, Table 21.2.2. Determines section classification from
    net tensile strain in the extreme tension reinforcement.

    Args:
        eps_t: Net tensile strain in extreme tension reinforcement.
        fy: Yield strength of reinforcement in MPa.
        Es: Modulus of elasticity in MPa (default 200000).
        transverse: Type of transverse reinforcement ('spiral' or 'other').

    Returns:
        Strength reduction factor (0.65 to 0.90).
    """
    eps_ty = fy / Es

    if eps_t <= eps_ty:
        return 0.75 if transverse == 'spiral' else 0.65
    elif eps_t >= eps_ty + 0.003:
        return 0.90
    else:
        if transverse == 'spiral':
            return 0.75 + 0.15 * (eps_t - eps_ty) / 0.003
        else:
            return 0.65 + 0.25 * (eps_t - eps_ty) / 0.003


def section_classification(
    eps_t: float,
    fy: float,
    Es: float = 200000.0,
) -> t.Literal['tension-controlled', 'transition', 'compression-controlled']:
    """Classify section per ACI 318-25, Table 21.2.2.

    Args:
        eps_t: Net tensile strain in extreme tension reinforcement.
        fy: Yield strength of reinforcement in MPa.
        Es: Modulus of elasticity in MPa (default 200000).

    Returns:
        Section classification string.
    """
    eps_ty = fy / Es
    if eps_t >= eps_ty + 0.003:
        return 'tension-controlled'
    elif eps_t <= eps_ty:
        return 'compression-controlled'
    else:
        return 'transition'
```

- [ ] **Step 4: Update `__init__.py` exports**

Add to `structuralcodes/codes/aci318_25/__init__.py`:

```python
from ._strength_reduction import (
    phi_bearing,
    phi_flexure,
    phi_shear,
    phi_torsion,
    section_classification,
)
```

And add to `__all__`:
```python
'phi_bearing',
'phi_flexure',
'phi_shear',
'phi_torsion',
'section_classification',
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_aci318_25/test_strength_reduction.py -v`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add structuralcodes/codes/aci318_25/_strength_reduction.py structuralcodes/codes/aci318_25/__init__.py tests/test_aci318_25/test_strength_reduction.py
git commit -m "feat(aci318_25): add strength reduction factors (Ch. 21)"
```

---

## Task 5: Flexure Functions

**Files:**
- Create: `structuralcodes/codes/aci318_25/_flexure.py`
- Create: `tests/test_aci318_25/test_flexure.py`
- Modify: `structuralcodes/codes/aci318_25/__init__.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_aci318_25/test_flexure.py
"""Tests for flexural strength functions of ACI 318-25, Ch. 22.2-22.3."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _flexure as fl

# Common parameters: 4000 psi concrete, Gr 60 steel, 12" strip
FC = 27.58  # MPa (4000 psi)
FY = 420.0  # MPa (60 ksi)
B = 305.0   # mm (12 in.)
D = 227.0   # mm (~8.94 in., for 10" slab with #5 bars, 3/4" cover)
H = 254.0   # mm (10 in.)
BETA1 = 0.85


class TestStressBlockDepthSR:
    def test_known_value(self):
        """As = 645 mm2 (2 #5 bars = 2 * 200 mm2 approx, use 645 for check).
        a = 645 * 420 / (0.85 * 27.58 * 305) = 37.9 mm"""
        As = 645.0
        a = fl.stress_block_depth_sr(As, FY, FC, B)
        expected = As * FY / (0.85 * FC * B)
        assert math.isclose(a, expected, rel_tol=1e-6)


class TestStressBlockDepthDR:
    def test_known_value(self):
        As = 800.0
        As_prime = 200.0
        a = fl.stress_block_depth_dr(As, As_prime, FY, FY, FC, B)
        expected = (As * FY - As_prime * FY) / (0.85 * FC * B)
        assert math.isclose(a, expected, rel_tol=1e-6)


class TestNeutralAxisDepth:
    def test_known_value(self):
        a = 37.9
        c = fl.neutral_axis_depth(a, BETA1)
        assert math.isclose(c, a / BETA1, rel_tol=1e-6)


class TestEpsTFromC:
    def test_typical_slab(self):
        """c = 44.6 mm, d = 227 mm, eps_cu = 0.003.
        eps_t = 0.003 * (227 - 44.6) / 44.6 = 0.01227"""
        c = 44.6
        eps_t = fl.eps_t_from_c(c, D)
        expected = 0.003 * (D - c) / c
        assert math.isclose(eps_t, expected, rel_tol=1e-6)


class TestEpsSPrime:
    def test_known_value(self):
        """c = 80 mm, d' = 40 mm.
        eps_s' = 0.003 * (80 - 40) / 80 = 0.0015"""
        eps = fl.eps_s_prime(80, 40)
        assert math.isclose(eps, 0.0015, rel_tol=1e-6)


class TestMnSinglyReinforced:
    def test_known_value(self):
        """As = 645 mm2, a = 37.9 mm.
        Mn = 645 * 420 * (227 - 37.9/2) = 56.36e6 N-mm"""
        As = 645.0
        Mn = fl.Mn_singly_reinforced(As, FY, FC, B, D)
        a = As * FY / (0.85 * FC * B)
        expected = As * FY * (D - a / 2)
        assert math.isclose(Mn, expected, rel_tol=1e-6)


class TestMnDoublyReinforced:
    def test_known_value(self):
        As = 800.0
        As_prime = 200.0
        d_prime = 40.0
        Mn = fl.Mn_doubly_reinforced(As, As_prime, FY, FY, FC, B, D, d_prime)
        a = (As * FY - As_prime * FY) / (0.85 * FC * B)
        expected = (As * FY - As_prime * FY) * (D - a / 2) + As_prime * FY * (D - d_prime)
        assert math.isclose(Mn, expected, rel_tol=1e-6)


class TestAsMinSlab:
    def test_gr60(self):
        """7.6.1.1 -> 24.4.3.2: 0.0018 * b * h for Gr 60."""
        As_min = fl.As_min_slab(FY, B, H)
        assert math.isclose(As_min, 0.0018 * B * H, rel_tol=1e-6)

    def test_gr40(self):
        """24.4.3.2: 0.0020 * b * h for Gr 40/50."""
        As_min = fl.As_min_slab(280.0, B, H)
        assert math.isclose(As_min, 0.0020 * B * H, rel_tol=1e-6)

    def test_gr80(self):
        """24.4.3.2: max(0.0014, 0.0018*60000/fy_psi) * b * h for Gr 80+."""
        fy_psi = 550 * 145.038  # ~79771 psi
        ratio = max(0.0014, 0.0018 * 60000 / fy_psi)
        As_min = fl.As_min_slab(550.0, B, H)
        assert math.isclose(As_min, ratio * B * H, rel_tol=1e-3)


class TestAsMinBeam:
    def test_known_value(self):
        """9.6.1.2: max(0.25*sqrt(f'c)/fy, 1.4/fy) * bw * d."""
        bw = 305.0
        As_min = fl.As_min_beam(FC, FY, bw, D)
        expected = max(0.25 * math.sqrt(FC) / FY, 1.4 / FY) * bw * D
        assert math.isclose(As_min, expected, rel_tol=1e-6)


class TestAsMaxCheck:
    def test_tension_controlled(self):
        assert fl.As_max_check(eps_t=0.010, fy=FY) is True

    def test_not_tension_controlled(self):
        assert fl.As_max_check(eps_t=0.002, fy=FY) is False


class TestAsRequired:
    def test_round_trip(self):
        """Compute As for a known Mu, then verify Mn matches."""
        # Use phi=0.9, Mu = 50e6 N-mm
        Mu = 50e6
        phi = 0.9
        As = fl.As_required(Mu, phi, FY, FC, B, D)
        # Verify: Mn = As * fy * (d - a/2)
        a = As * FY / (0.85 * FC * B)
        Mn = As * FY * (D - a / 2)
        assert math.isclose(phi * Mn, Mu, rel_tol=1e-4)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_aci318_25/test_flexure.py -v`

Expected: ERRORS — module does not exist.

- [ ] **Step 3: Implement flexure functions**

```python
# structuralcodes/codes/aci318_25/_flexure.py
"""Flexural strength functions according to ACI 318-25, Ch. 22.2-22.3."""

from __future__ import annotations

import math


def stress_block_depth_sr(
    As: float, fy: float, fc: float, b: float,
) -> float:
    """Stress block depth a for singly-reinforced rectangular section.

    From equilibrium: a = As * fy / (0.85 * f'c * b).

    Args:
        As: Area of tension reinforcement (mm2).
        fy: Yield strength of reinforcement (MPa).
        fc: Specified compressive strength f'c (MPa).
        b: Width of compression face (mm).

    Returns:
        Stress block depth a (mm).
    """
    return As * fy / (0.85 * fc * b)


def stress_block_depth_dr(
    As: float, As_prime: float,
    fy: float, fy_prime: float,
    fc: float, b: float,
) -> float:
    """Stress block depth a for doubly-reinforced rectangular section.

    From equilibrium: a = (As*fy - As'*fy') / (0.85 * f'c * b).

    Args:
        As: Area of tension reinforcement (mm2).
        As_prime: Area of compression reinforcement (mm2).
        fy: Yield strength of tension reinforcement (MPa).
        fy_prime: Stress in compression reinforcement (MPa).
        fc: Specified compressive strength f'c (MPa).
        b: Width of compression face (mm).

    Returns:
        Stress block depth a (mm).
    """
    return (As * fy - As_prime * fy_prime) / (0.85 * fc * b)


def neutral_axis_depth(a: float, beta1: float) -> float:
    """Neutral axis depth from stress block depth.

    c = a / beta1.

    Args:
        a: Stress block depth (mm).
        beta1: Stress block depth factor.

    Returns:
        Neutral axis depth c (mm).
    """
    return a / beta1


def eps_t_from_c(
    c: float, dt: float, eps_cu: float = 0.003,
) -> float:
    """Net tensile strain in extreme tension reinforcement.

    ACI 318-25, Fig. R21.2.2a.
    eps_t = eps_cu * (dt - c) / c

    Args:
        c: Neutral axis depth (mm).
        dt: Distance from extreme compression fiber to extreme
            tension reinforcement (mm).
        eps_cu: Ultimate concrete strain (default 0.003).

    Returns:
        Net tensile strain (dimensionless).
    """
    return eps_cu * (dt - c) / c


def eps_s_prime(
    c: float, d_prime: float, eps_cu: float = 0.003,
) -> float:
    """Strain in compression reinforcement.

    eps_s' = eps_cu * (c - d') / c

    Args:
        c: Neutral axis depth (mm).
        d_prime: Distance from extreme compression fiber to
            compression reinforcement (mm).
        eps_cu: Ultimate concrete strain (default 0.003).

    Returns:
        Compression steel strain (dimensionless).
    """
    return eps_cu * (c - d_prime) / c


def Mn_singly_reinforced(
    As: float, fy: float, fc: float, b: float, d: float,
) -> float:
    """Nominal flexural strength of singly-reinforced rectangular section.

    ACI 318-25, Sec. 22.3.
    Mn = As * fy * (d - a/2)

    Args:
        As: Area of tension reinforcement (mm2).
        fy: Yield strength of reinforcement (MPa).
        fc: Specified compressive strength f'c (MPa).
        b: Width of compression face (mm).
        d: Effective depth (mm).

    Returns:
        Nominal moment strength Mn (N-mm).
    """
    a = stress_block_depth_sr(As, fy, fc, b)
    return As * fy * (d - a / 2)


def Mn_doubly_reinforced(
    As: float, As_prime: float,
    fy: float, fy_prime: float,
    fc: float, b: float, d: float, d_prime: float,
) -> float:
    """Nominal flexural strength of doubly-reinforced rectangular section.

    Mn = (As*fy - As'*fy') * (d - a/2) + As'*fy' * (d - d')

    The caller must verify compression steel yields via eps_s_prime().
    If it hasn't yielded, fy_prime should be replaced with Es * eps_s'.

    Args:
        As: Area of tension reinforcement (mm2).
        As_prime: Area of compression reinforcement (mm2).
        fy: Yield strength of tension reinforcement (MPa).
        fy_prime: Stress in compression reinforcement (MPa).
        fc: Specified compressive strength f'c (MPa).
        b: Width of compression face (mm).
        d: Effective depth to tension reinforcement (mm).
        d_prime: Depth to compression reinforcement (mm).

    Returns:
        Nominal moment strength Mn (N-mm).
    """
    a = stress_block_depth_dr(As, As_prime, fy, fy_prime, fc, b)
    return (As * fy - As_prime * fy_prime) * (d - a / 2) + (
        As_prime * fy_prime * (d - d_prime)
    )


def As_min_slab(fy: float, b: float, h: float) -> float:
    """Minimum flexural reinforcement for one-way slabs.

    ACI 318-25, Sec. 7.6.1.1 -> 24.4.3.2.
    Same as shrinkage and temperature reinforcement.

    Args:
        fy: Yield strength of reinforcement (MPa).
        b: Width of slab strip (mm).
        h: Overall slab thickness (mm).

    Returns:
        Minimum reinforcement area (mm2).
    """
    fy_psi = fy * 145.038
    if fy_psi <= 50000:
        ratio = 0.0020
    elif fy_psi <= 60000:
        ratio = 0.0018
    else:
        ratio = max(0.0014, 0.0018 * 60000 / fy_psi)
    return ratio * b * h


def As_min_beam(
    fc: float, fy: float, bw: float, d: float,
) -> float:
    """Minimum flexural reinforcement for beams.

    ACI 318-25, Sec. 9.6.1.2.
    As_min = max(0.25*sqrt(f'c)/fy, 1.4/fy) * bw * d

    Args:
        fc: Specified compressive strength f'c (MPa).
        fy: Yield strength of reinforcement (MPa).
        bw: Web width (mm).
        d: Effective depth (mm).

    Returns:
        Minimum reinforcement area (mm2).
    """
    return max(0.25 * math.sqrt(fc) / fy, 1.4 / fy) * bw * d


def As_max_check(
    eps_t: float, fy: float, Es: float = 200000.0,
) -> bool:
    """Check that section is tension-controlled.

    Required for slabs (7.3.3.1) and beams (9.3.3.1).
    Tension-controlled: eps_t >= eps_ty + 0.003.

    Args:
        eps_t: Net tensile strain in extreme tension reinforcement.
        fy: Yield strength of reinforcement (MPa).
        Es: Modulus of elasticity (MPa).

    Returns:
        True if tension-controlled.
    """
    eps_ty = fy / Es
    return eps_t >= eps_ty + 0.003


def As_required(
    Mu: float, phi: float, fy: float, fc: float, b: float, d: float,
) -> float:
    """Required tension reinforcement area for singly-reinforced section.

    Solves: Mu = phi * As * fy * (d - As*fy / (1.7*f'c*b))
    via the quadratic formula.

    Args:
        Mu: Factored moment (N-mm).
        phi: Strength reduction factor.
        fy: Yield strength of reinforcement (MPa).
        fc: Specified compressive strength f'c (MPa).
        b: Width of compression face (mm).
        d: Effective depth (mm).

    Returns:
        Required reinforcement area (mm2).

    Raises:
        ValueError: If no real solution exists (section too small).
    """
    Mn_req = Mu / phi
    # Mn = rho*fy*b*d^2*(1 - 0.5*rho*fy/(0.85*fc))
    # Rearranging: 0 = As^2 * fy/(1.7*fc*b) - As * fy * d + Mn_req
    a_coeff = fy / (1.7 * fc * b)
    b_coeff = -fy * d
    c_coeff = Mn_req
    discriminant = b_coeff**2 - 4 * a_coeff * c_coeff
    if discriminant < 0:
        raise ValueError(
            'No real solution: section is too small for the required moment.'
        )
    return (-b_coeff - math.sqrt(discriminant)) / (2 * a_coeff)
```

- [ ] **Step 4: Update `__init__.py` exports**

Add to `structuralcodes/codes/aci318_25/__init__.py`:

```python
from ._flexure import (
    As_max_check,
    As_min_beam,
    As_min_slab,
    As_required,
    Mn_doubly_reinforced,
    Mn_singly_reinforced,
    eps_s_prime,
    eps_t_from_c,
    neutral_axis_depth,
    stress_block_depth_dr,
    stress_block_depth_sr,
)
```

And add all names to `__all__`.

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_aci318_25/test_flexure.py -v`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add structuralcodes/codes/aci318_25/_flexure.py structuralcodes/codes/aci318_25/__init__.py tests/test_aci318_25/test_flexure.py
git commit -m "feat(aci318_25): add flexural strength functions (Ch. 22.2-22.3)"
```

---

## Task 6: Shear Functions

**Files:**
- Create: `structuralcodes/codes/aci318_25/_shear.py`
- Create: `tests/test_aci318_25/test_shear.py`
- Modify: `structuralcodes/codes/aci318_25/__init__.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_aci318_25/test_shear.py
"""Tests for one-way shear strength functions of ACI 318-25, Ch. 22.5."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _shear as sh

FC = 27.58   # MPa (4000 psi)
BW = 305.0   # mm (12 in.)
D = 227.0    # mm
RHO_W = 0.009  # typical slab reinforcement ratio


class TestLambdaS:
    def test_small_member(self):
        """d = 227 mm = 8.94 in. lambda_s = 2/(1+8.94/10) = 1.056 -> capped at 1.0."""
        assert sh.lambda_s(227) == 1.0

    def test_large_member(self):
        """d = 900 mm = 35.4 in. lambda_s = 2/(1+35.4/10) = 0.440."""
        result = sh.lambda_s(900)
        d_in = 900 / 25.4
        expected = min(2 / (1 + d_in / 10), 1.0)
        assert math.isclose(result, expected, rel_tol=1e-3)


class TestVcDetailed:
    def test_without_min_reinforcement(self):
        """Table 22.5.5.1(c): Vc = 8*lambda_s*lambda*(rho_w)^(1/3)*sqrt(f'c)*bw*d."""
        Vc = sh.Vc_detailed(FC, BW, D, RHO_W, Av_provided=0, Av_min=100)
        ls = min(2 / (1 + (D / 25.4) / 10), 1.0)
        expected = 8 * ls * 1.0 * RHO_W**(1/3) * math.sqrt(FC) * BW * D
        # Apply upper bound
        upper = 5 * 1.0 * math.sqrt(FC) * BW * D
        expected = min(expected, upper)
        # Apply lower bound
        lower = 1.0 * math.sqrt(FC) * BW * D
        expected = max(expected, lower)
        assert math.isclose(Vc, expected, rel_tol=1e-3)

    def test_with_min_reinforcement(self):
        """Table 22.5.5.1(b): Vc = 8*lambda*(rho_w)^(1/3)*sqrt(f'c)*bw*d."""
        Vc = sh.Vc_detailed(FC, BW, D, RHO_W, Av_provided=200, Av_min=100)
        expected = 8 * 1.0 * RHO_W**(1/3) * math.sqrt(FC) * BW * D
        upper = 5 * 1.0 * math.sqrt(FC) * BW * D
        expected = min(expected, upper)
        lower = 1.0 * math.sqrt(FC) * BW * D
        expected = max(expected, lower)
        assert math.isclose(Vc, expected, rel_tol=1e-3)

    def test_vc_not_negative(self):
        """Vc with large axial tension should be >= 0."""
        Vc = sh.Vc_detailed(FC, BW, D, RHO_W, Nu=-500000, Ag=BW * 254)
        assert Vc >= 0


class TestVcSimplified:
    def test_no_axial(self):
        """Table 22.5.5.1(a): Vc = 2*lambda*sqrt(f'c)*bw*d."""
        Vc = sh.Vc_simplified(FC, BW, D)
        expected = 2 * 1.0 * math.sqrt(FC) * BW * D
        assert math.isclose(Vc, expected, rel_tol=1e-6)


class TestVs:
    def test_known_value(self):
        """Av=142 mm2 (2 legs #3), fyt=420 MPa, d=227 mm, s=150 mm."""
        result = sh.Vs(142, 420, 227, 150)
        expected = 142 * 420 * 227 / 150
        assert math.isclose(result, expected, rel_tol=1e-6)


class TestVn:
    def test_sum(self):
        assert sh.Vn(50000, 30000) == 80000


class TestCheckCrossSection:
    def test_passes(self):
        Vc = 50000
        assert sh.check_cross_section(30000, 0.75, Vc, FC, BW, D) is True

    def test_fails(self):
        Vc = 50000
        huge_Vu = 1e7
        assert sh.check_cross_section(huge_Vu, 0.75, Vc, FC, BW, D) is False


class TestAvMinPerS:
    def test_known_value(self):
        """max(0.062*sqrt(f'c), 0.35) * bw / fyt."""
        fyt = 420.0
        result = sh.Av_min_per_s(FC, BW, fyt)
        expected = max(0.062 * math.sqrt(FC), 0.35) * BW / fyt
        assert math.isclose(result, expected, rel_tol=1e-6)


class TestShearReinforcementRequired:
    def test_not_required(self):
        assert sh.shear_reinforcement_required(30000, 50000) is False

    def test_required(self):
        assert sh.shear_reinforcement_required(60000, 50000) is True


class TestMaxStirrupSpacing:
    def test_low_vs(self):
        """Vs <= 4*sqrt(f'c)*bw*d -> s_max = min(d/2, 600)."""
        Vs = 1000  # very low
        s_max = sh.max_stirrup_spacing(D, Vs, FC, BW)
        assert math.isclose(s_max, min(D / 2, 600), rel_tol=1e-6)

    def test_high_vs(self):
        """Vs > 4*sqrt(f'c)*bw*d -> s_max = min(d/4, 300)."""
        Vs = 4 * math.sqrt(FC) * BW * D + 1  # just above threshold
        s_max = sh.max_stirrup_spacing(D, Vs, FC, BW)
        assert math.isclose(s_max, min(D / 4, 300), rel_tol=1e-6)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_aci318_25/test_shear.py -v`

Expected: ERRORS — module does not exist.

- [ ] **Step 3: Implement shear functions**

```python
# structuralcodes/codes/aci318_25/_shear.py
"""One-way shear strength functions according to ACI 318-25, Ch. 22.5."""

from __future__ import annotations

import math


def lambda_s(d: float) -> float:
    """Size effect modification factor.

    ACI 318-25, Eq. 22.5.5.1.3.
    lambda_s = 2 / (1 + d_in/10) <= 1.0

    where d_in is in inches. This function accepts d in mm.

    Args:
        d: Effective depth in mm.

    Returns:
        Size effect factor (dimensionless), <= 1.0.
    """
    d_in = d / 25.4
    return min(2.0 / (1.0 + d_in / 10.0), 1.0)


def Vc_detailed(
    fc: float,
    bw: float,
    d: float,
    rho_w: float,
    Nu: float = 0.0,
    Ag: float = 0.0,
    lambda_concrete: float = 1.0,
    Av_provided: float = 0.0,
    Av_min: float = 0.0,
) -> float:
    """Concrete shear strength for nonprestressed members.

    ACI 318-25, Table 22.5.5.1.

    If Av >= Av_min (b): Vc = [8*lambda*(rho_w)^(1/3)*sqrt(f'c) + Nu/(6*Ag)] * bw * d
    If Av < Av_min  (c): Vc = [8*lambda_s*lambda*(rho_w)^(1/3)*sqrt(f'c) + Nu/(6*Ag)] * bw * d

    Limits:
        Vc <= 5*lambda*sqrt(f'c)*bw*d  (22.5.5.1.1)
        Vc >= lambda*sqrt(f'c)*bw*d    (22.5.5.1.1, unless net axial tension)
        Nu/(6*Ag) <= 0.05*f'c          (22.5.5.1.2)
        sqrt(f'c) <= 8.3 MPa           (22.5.3.1, ~100 psi)

    Args:
        fc: Specified compressive strength f'c (MPa).
        bw: Web width (mm).
        d: Effective depth (mm).
        rho_w: Longitudinal reinforcement ratio As/(bw*d).
        Nu: Axial force (N), positive for compression, negative for tension.
        Ag: Gross area (mm2).
        lambda_concrete: Lightweight modification factor (default 1.0).
        Av_provided: Provided shear reinforcement area per spacing (mm2/mm).
        Av_min: Required minimum shear reinforcement per spacing (mm2/mm).

    Returns:
        Concrete shear strength Vc (N).
    """
    sqrt_fc = min(math.sqrt(fc), 8.3)
    lam = lambda_concrete

    # Axial load term (22.5.5.1.2)
    axial_term = 0.0
    if Ag > 0:
        axial_term = min(Nu / (6.0 * Ag), 0.05 * fc)

    # Reinforcement ratio term
    rho_term = 8.0 * lam * rho_w ** (1.0 / 3.0) * sqrt_fc

    if Av_provided >= Av_min:
        vc = rho_term + axial_term
    else:
        ls = lambda_s(d)
        vc = 8.0 * ls * lam * rho_w ** (1.0 / 3.0) * sqrt_fc + axial_term

    Vc = vc * bw * d

    # Upper bound (22.5.5.1.1)
    Vc_max = 5.0 * lam * sqrt_fc * bw * d
    Vc = min(Vc, Vc_max)

    # Lower bound (22.5.5.1.1) — does not apply for net axial tension
    if Nu >= 0:
        Vc_min = lam * sqrt_fc * bw * d
        Vc = max(Vc, Vc_min)

    # Vc shall not be less than zero (Table 22.5.5.1, Note 2)
    return max(Vc, 0.0)


def Vc_simplified(
    fc: float,
    bw: float,
    d: float,
    Nu: float = 0.0,
    Ag: float = 0.0,
    lambda_concrete: float = 1.0,
) -> float:
    """Simplified concrete shear strength.

    ACI 318-25, Table 22.5.5.1(a).
    Vc = [2*lambda*sqrt(f'c) + Nu/(6*Ag)] * bw * d

    Only valid when Av >= Av_min.

    Args:
        fc: Specified compressive strength f'c (MPa).
        bw: Web width (mm).
        d: Effective depth (mm).
        Nu: Axial force (N), positive for compression.
        Ag: Gross area (mm2).
        lambda_concrete: Lightweight modification factor (default 1.0).

    Returns:
        Concrete shear strength Vc (N).
    """
    sqrt_fc = min(math.sqrt(fc), 8.3)
    axial_term = 0.0
    if Ag > 0:
        axial_term = min(Nu / (6.0 * Ag), 0.05 * fc)
    Vc = (2.0 * lambda_concrete * sqrt_fc + axial_term) * bw * d
    return max(Vc, 0.0)


def Vs(Av: float, fyt: float, d: float, s: float) -> float:
    """Shear strength provided by transverse reinforcement.

    ACI 318-25, Eq. 22.5.8.5.3.

    Args:
        Av: Area of shear reinforcement within spacing s (mm2).
        fyt: Yield strength of transverse reinforcement (MPa).
        d: Effective depth (mm).
        s: Spacing of transverse reinforcement (mm).

    Returns:
        Shear strength Vs (N).
    """
    return Av * fyt * d / s


def Vn(Vc: float, Vs: float) -> float:
    """Nominal shear strength.

    Args:
        Vc: Concrete contribution (N).
        Vs: Steel contribution (N).

    Returns:
        Nominal shear strength (N).
    """
    return Vc + Vs


def check_cross_section(
    Vu: float, phi: float, Vc: float, fc: float, bw: float, d: float,
) -> bool:
    """Check cross-section dimensions.

    ACI 318-25, Eq. 22.5.1.2.
    Vu <= phi * (Vc + 8*sqrt(f'c)*bw*d)

    Returns:
        True if dimensions are adequate.
    """
    sqrt_fc = min(math.sqrt(fc), 8.3)
    return Vu <= phi * (Vc + 8.0 * sqrt_fc * bw * d)


def Av_min_per_s(fc: float, bw: float, fyt: float) -> float:
    """Minimum shear reinforcement area per unit spacing.

    ACI 318-25, Sec. 9.6.3.4.
    Av_min/s = max(0.062*sqrt(f'c), 0.35) * bw / fyt

    Args:
        fc: Specified compressive strength f'c (MPa).
        bw: Web width (mm).
        fyt: Yield strength of transverse reinforcement (MPa).

    Returns:
        Minimum Av/s (mm2/mm).
    """
    return max(0.062 * math.sqrt(fc), 0.35) * bw / fyt


def shear_reinforcement_required(Vu: float, phi_Vc: float) -> bool:
    """Whether shear reinforcement is required.

    ACI 318-25, Sec. 7.6.3.1: required when Vu > phi*Vc.

    Returns:
        True if shear reinforcement is required.
    """
    return Vu > phi_Vc


def max_stirrup_spacing(
    d: float, Vs: float, fc: float, bw: float,
) -> float:
    """Maximum stirrup spacing.

    ACI 318-25, Sec. 9.7.6.2.2.
    If Vs <= 4*sqrt(f'c)*bw*d: s_max = min(d/2, 600 mm)
    If Vs >  4*sqrt(f'c)*bw*d: s_max = min(d/4, 300 mm)

    Returns:
        Maximum spacing (mm).
    """
    threshold = 4.0 * math.sqrt(fc) * bw * d
    if Vs <= threshold:
        return min(d / 2, 600.0)
    return min(d / 4, 300.0)
```

- [ ] **Step 4: Update `__init__.py` exports**

Add to `structuralcodes/codes/aci318_25/__init__.py`:

```python
from ._shear import (
    Av_min_per_s,
    Vc_detailed,
    Vc_simplified,
    Vn,
    Vs,
    check_cross_section,
    lambda_s,
    max_stirrup_spacing,
    shear_reinforcement_required,
)
```

And add all names to `__all__`.

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_aci318_25/test_shear.py -v`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add structuralcodes/codes/aci318_25/_shear.py structuralcodes/codes/aci318_25/__init__.py tests/test_aci318_25/test_shear.py
git commit -m "feat(aci318_25): add one-way shear strength functions (Ch. 22.5)"
```

---

## Task 7: One-Way Slab Module

**Files:**
- Create: `structuralcodes/codes/aci318_25/_one_way_slab.py`
- Create: `tests/test_aci318_25/test_one_way_slab.py`
- Modify: `structuralcodes/codes/aci318_25/__init__.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_aci318_25/test_one_way_slab.py
"""Tests for one-way slab rules of ACI 318-25, Chapter 7."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _one_way_slab as ows


class TestMinThickness:
    def test_simply_supported_gr60(self):
        """L/20 for simply supported, fy=420 MPa (Gr 60)."""
        span = 6096.0  # 20 ft in mm
        h = ows.min_thickness(span, 'simply_supported')
        assert math.isclose(h, span / 20, rel_tol=1e-6)

    def test_one_end_continuous_gr60(self):
        """L/24 for one end continuous."""
        span = 6096.0
        h = ows.min_thickness(span, 'one_end_continuous')
        assert math.isclose(h, span / 24, rel_tol=1e-6)

    def test_both_ends_continuous_gr60(self):
        """L/28 for both ends continuous."""
        span = 6096.0
        h = ows.min_thickness(span, 'both_ends_continuous')
        assert math.isclose(h, span / 28, rel_tol=1e-6)

    def test_cantilever_gr60(self):
        """L/10 for cantilever."""
        span = 3048.0  # 10 ft
        h = ows.min_thickness(span, 'cantilever')
        assert math.isclose(h, span / 10, rel_tol=1e-6)

    def test_fy_adjustment(self):
        """7.3.1.1.1: multiply by (0.4 + fy/100000) for fy != 60 ksi."""
        span = 6096.0
        fy_80 = 550.0  # ~80 ksi
        h_60 = ows.min_thickness(span, 'simply_supported', fy=420.0)
        h_80 = ows.min_thickness(span, 'simply_supported', fy=550.0)
        fy_psi = 550 * 145.038
        factor = 0.4 + fy_psi / 100000
        assert math.isclose(h_80, h_60 * factor, rel_tol=1e-3)

    def test_invalid_support_raises(self):
        with pytest.raises(ValueError):
            ows.min_thickness(6096, 'invalid')


class TestAsShrinkageTemperature:
    def test_gr60(self):
        """0.0018 * b * h for Gr 60."""
        b, h = 305.0, 254.0
        As = ows.As_shrinkage_temperature(420, b, h)
        assert math.isclose(As, 0.0018 * b * h, rel_tol=1e-6)


class TestMaxBarSpacingFlexure:
    def test_thin_slab(self):
        """min(3*h, 450). For h=150, 3*150=450."""
        assert math.isclose(ows.max_bar_spacing_flexure(150), 450, rel_tol=1e-6)

    def test_thick_slab(self):
        """For h=200, 3*200=600 > 450, so 450."""
        assert math.isclose(ows.max_bar_spacing_flexure(200), 450, rel_tol=1e-6)


class TestMaxBarSpacingShrinkage:
    def test_value(self):
        """min(5*h, 450). For h=100, 5*100=500 > 450, so 450."""
        assert math.isclose(ows.max_bar_spacing_shrinkage(100), 450, rel_tol=1e-6)

    def test_thin_slab(self):
        """For h=80, 5*80=400 < 450, so 400."""
        assert math.isclose(ows.max_bar_spacing_shrinkage(80), 400, rel_tol=1e-6)


class TestShearCriticalSectionOffset:
    def test_value(self):
        assert ows.shear_critical_section_offset(227) == 227
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_aci318_25/test_one_way_slab.py -v`

Expected: ERRORS — module does not exist.

- [ ] **Step 3: Implement one-way slab functions**

```python
# structuralcodes/codes/aci318_25/_one_way_slab.py
"""One-way slab design rules according to ACI 318-25, Chapter 7."""

from __future__ import annotations

import typing as t

THICKNESS_RATIOS = {
    'simply_supported': 20,
    'one_end_continuous': 24,
    'both_ends_continuous': 28,
    'cantilever': 10,
}


def min_thickness(
    span: float,
    support_condition: t.Literal[
        'simply_supported', 'one_end_continuous',
        'both_ends_continuous', 'cantilever',
    ],
    fy: float = 420.0,
    lightweight: bool = False,
    wc: float = 2320.0,
) -> float:
    """Minimum slab thickness to satisfy deflection without calculation.

    ACI 318-25, Table 7.3.1.1.

    Base ratios: L/20 (simply supported), L/24 (one end continuous),
    L/28 (both ends continuous), L/10 (cantilever).

    Adjustments:
        7.3.1.1.1: For fy != 60 ksi, multiply by (0.4 + fy/100000) [fy in psi].
        7.3.1.1.2: For lightweight, multiply by max(1.65 - 0.005*wc, 1.09) [wc in pcf].

    Args:
        span: Clear span length (mm).
        support_condition: End condition string.
        fy: Yield strength (MPa, default 420 for Gr 60).
        lightweight: Whether lightweight concrete.
        wc: Unit weight (kg/m3), only used if lightweight=True.

    Returns:
        Minimum thickness h (mm).

    Raises:
        ValueError: If support_condition is not recognized.
    """
    ratio = THICKNESS_RATIOS.get(support_condition)
    if ratio is None:
        raise ValueError(
            f'Unknown support condition: {support_condition}. '
            f'Valid options: {list(THICKNESS_RATIOS.keys())}'
        )

    h = span / ratio

    # Adjustment for fy (7.3.1.1.1)
    fy_psi = fy * 145.038
    if not (59000 < fy_psi < 61000):
        h *= (0.4 + fy_psi / 100000)

    # Adjustment for lightweight (7.3.1.1.2)
    if lightweight:
        wc_pcf = wc / 16.0185
        h *= max(1.65 - 0.005 * wc_pcf, 1.09)

    return h


def As_shrinkage_temperature(fy: float, b: float, h: float) -> float:
    """Shrinkage and temperature reinforcement area.

    ACI 318-25, Sec. 24.4.3.2.
    Also the minimum flexural reinforcement for slabs (7.6.1.1).

    Args:
        fy: Yield strength of reinforcement (MPa).
        b: Width of slab strip (mm).
        h: Overall slab thickness (mm).

    Returns:
        Required reinforcement area (mm2).
    """
    fy_psi = fy * 145.038
    if fy_psi <= 50000:
        ratio = 0.0020
    elif fy_psi <= 60000:
        ratio = 0.0018
    else:
        ratio = max(0.0014, 0.0018 * 60000 / fy_psi)
    return ratio * b * h


def max_bar_spacing_flexure(h: float) -> float:
    """Maximum spacing of flexural reinforcement.

    ACI 318-25, Sec. 7.7.2.3.
    s_max = min(3*h, 450 mm)

    Args:
        h: Overall slab thickness (mm).

    Returns:
        Maximum spacing (mm).
    """
    return min(3.0 * h, 450.0)


def max_bar_spacing_shrinkage(h: float) -> float:
    """Maximum spacing of shrinkage/temperature reinforcement.

    ACI 318-25, Sec. 7.7.6.2.1.
    s_max = min(5*h, 450 mm)

    Args:
        h: Overall slab thickness (mm).

    Returns:
        Maximum spacing (mm).
    """
    return min(5.0 * h, 450.0)


def shear_critical_section_offset(d: float) -> float:
    """Distance from face of support to critical section for shear.

    ACI 318-25, Sec. 7.4.3.2.
    For nonprestressed slabs: d from face of support.

    Args:
        d: Effective depth (mm).

    Returns:
        Offset distance (mm).
    """
    return d
```

- [ ] **Step 4: Update `__init__.py` exports**

Add to `structuralcodes/codes/aci318_25/__init__.py`:

```python
from ._one_way_slab import (
    As_shrinkage_temperature,
    max_bar_spacing_flexure,
    max_bar_spacing_shrinkage,
    min_thickness,
    shear_critical_section_offset,
)
```

And add all names to `__all__`.

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_aci318_25/test_one_way_slab.py -v`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add structuralcodes/codes/aci318_25/_one_way_slab.py structuralcodes/codes/aci318_25/__init__.py tests/test_aci318_25/test_one_way_slab.py
git commit -m "feat(aci318_25): add one-way slab design rules (Ch. 7)"
```

---

## Task 8: ConcreteACI318_25 Material Class

**Files:**
- Create: `structuralcodes/materials/concrete/_concreteACI318_25.py`
- Create: `tests/test_aci318_25/test_concrete_aci318_25.py`
- Modify: `structuralcodes/materials/concrete/__init__.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_aci318_25/test_concrete_aci318_25.py
"""Tests for the ConcreteACI318_25 material class."""

import math

import pytest

import structuralcodes
from structuralcodes.materials.concrete._concreteACI318_25 import ConcreteACI318_25


@pytest.fixture(autouse=True)
def _reset_design_code():
    yield
    structuralcodes.set_design_code(None)


class TestConstruction:
    def test_basic(self):
        c = ConcreteACI318_25(fck=27.58)
        assert c.fck == 27.58

    def test_fc_alias(self):
        c = ConcreteACI318_25(fck=27.58)
        assert c.fc == 27.58

    def test_default_name(self):
        c = ConcreteACI318_25(fck=28)
        assert c.name == 'C28'

    def test_custom_name(self):
        c = ConcreteACI318_25(fck=28, name='4000psi')
        assert c.name == '4000psi'


class TestProperties:
    def test_gamma_c(self):
        c = ConcreteACI318_25(fck=27.58)
        assert c.gamma_c == 1.0

    def test_fcd(self):
        c = ConcreteACI318_25(fck=27.58)
        assert math.isclose(c.fcd(), 0.85 * 27.58, rel_tol=1e-6)

    def test_Ec(self):
        c = ConcreteACI318_25(fck=27.58)
        expected = 2320**1.5 * 0.043 * math.sqrt(27.58)
        assert math.isclose(c.Ec, expected, rel_tol=5e-3)

    def test_Ec_override(self):
        c = ConcreteACI318_25(fck=27.58, Ec=25000)
        assert c.Ec == 25000

    def test_fr(self):
        c = ConcreteACI318_25(fck=27.58)
        assert math.isclose(c.fr, 0.62 * math.sqrt(27.58), rel_tol=1e-6)

    def test_beta1(self):
        c = ConcreteACI318_25(fck=27.58)
        assert c.beta1 == 0.85

    def test_alpha1(self):
        c = ConcreteACI318_25(fck=27.58)
        assert c.alpha1 == 0.85

    def test_eps_cu(self):
        c = ConcreteACI318_25(fck=27.58)
        assert c.eps_cu == 0.003


class TestConstitutiveLaws:
    def test_elastic(self):
        c = ConcreteACI318_25(fck=27.58, constitutive_law='elastic')
        assert c.constitutive_law is not None

    def test_parabolarectangle(self):
        c = ConcreteACI318_25(fck=27.58, constitutive_law='parabolarectangle')
        assert c.constitutive_law is not None
        # Check ultimate strain is 0.003
        eps_min, eps_max = c.constitutive_law.get_ultimate_strain()
        assert math.isclose(abs(eps_min), 0.003, rel_tol=1e-6)

    def test_bilinearcompression(self):
        c = ConcreteACI318_25(fck=27.58, constitutive_law='bilinearcompression')
        assert c.constitutive_law is not None


class TestFactory:
    def test_create_via_factory(self):
        from structuralcodes.materials.concrete import create_concrete
        c = create_concrete(fck=27.58, design_code='aci318_25')
        assert isinstance(c, ConcreteACI318_25)

    def test_create_via_global_code(self):
        from structuralcodes.materials.concrete import create_concrete
        structuralcodes.set_design_code('aci318_25')
        c = create_concrete(fck=27.58)
        assert isinstance(c, ConcreteACI318_25)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_aci318_25/test_concrete_aci318_25.py -v`

Expected: ERRORS — module does not exist.

- [ ] **Step 3: Implement ConcreteACI318_25**

```python
# structuralcodes/materials/concrete/_concreteACI318_25.py
"""Concrete material class for ACI 318-25."""

import typing as t

from structuralcodes.codes import aci318_25

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._concrete import Concrete


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

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400,
        gamma_c: t.Optional[float] = None,
        constitutive_law: t.Optional[
            t.Union[
                t.Literal[
                    'elastic',
                    'parabolarectangle',
                    'bilinearcompression',
                ],
                ConstitutiveLaw,
            ]
        ] = 'parabolarectangle',
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[bool] = None,
        Ec: t.Optional[float] = None,
        fr: t.Optional[float] = None,
        wc: float = 2320.0,
        lambda_s: float = 1.0,
        **kwargs,
    ) -> None:
        """Initialize ACI 318-25 concrete.

        Args:
            fck: Specified compressive strength f'c in MPa.
            name: Descriptive name (default: 'C{fck}').
            density: Density in kg/m3 (default 2400).
            gamma_c: Partial factor (default 1.0 for ACI).
            constitutive_law: ConstitutiveLaw or string name.
            initial_strain: Initial strain of the material.
            initial_stress: Initial stress of the material.
            strain_compatibility: Whether material deforms with geometry.
            Ec: Override for modulus of elasticity (MPa).
            fr: Override for modulus of rupture (MPa).
            wc: Unit weight for Ec calculation (kg/m3, default 2320).
            lambda_s: Lightweight modification factor (default 1.0).
        """
        del kwargs
        if name is None:
            name = f'C{round(fck):d}'
        super().__init__(
            fck=fck,
            name=name,
            density=density,
            existing=False,
            gamma_c=gamma_c,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            strain_compatibility=strain_compatibility,
        )
        self._Ec = Ec
        self._fr = fr
        self._wc = wc
        self._lambda_s = lambda_s

        self._constitutive_law = (
            constitutive_law
            if isinstance(constitutive_law, ConstitutiveLaw)
            else create_constitutive_law(
                constitutive_law_name=constitutive_law, material=self
            )
        )
        if 'concrete' not in self._constitutive_law.__materials__:
            raise ValueError(
                'The provided constitutive law is not valid for concrete.'
            )
        self._apply_initial_strain()

    @property
    def fc(self) -> float:
        """Specified compressive strength f'c in MPa (alias for fck)."""
        return self._fck

    @property
    def gamma_c(self) -> float:
        """Partial factor for concrete.

        Returns 1.0 for ACI 318. ACI does not reduce material strengths;
        safety is applied via phi factors at member capacity level.
        """
        return self._gamma_c or 1.0

    @property
    def alpha1(self) -> float:
        """Stress block intensity factor (Sec. 22.2.2.4.1)."""
        return aci318_25.alpha1()

    def fcd(self) -> float:
        """Design compressive strength.

        Returns alpha1 * f'c / gamma_c = 0.85 * f'c for ACI.
        This is the stress intensity for the Whitney stress block,
        not a gamma-reduced design strength.
        """
        return self.alpha1 * self.fc / self.gamma_c

    @property
    def Ec(self) -> float:
        """Modulus of elasticity (Table 19.2.2.1) in MPa."""
        if self._Ec is not None:
            return self._Ec
        return aci318_25.Ec(self.fc, wc=self._wc)

    @property
    def fr(self) -> float:
        """Modulus of rupture (Eq. 19.2.3.1) in MPa."""
        if self._fr is not None:
            return self._fr
        return aci318_25.fr(self.fc, lambda_s=self._lambda_s)

    @property
    def fct(self) -> float:
        """Splitting tensile strength (Sec. 19.2.4.3) in MPa."""
        return aci318_25.fct(self.fc, lambda_s=self._lambda_s)

    @property
    def beta1(self) -> float:
        """Whitney stress block depth factor (Table 22.2.2.4.3)."""
        return aci318_25.beta1(self.fc)

    @property
    def eps_cu(self) -> float:
        """Ultimate concrete strain (Sec. 22.2.2.1)."""
        return aci318_25.eps_cu()

    def __elastic__(self) -> dict:
        """Returns kwargs for creating an elastic constitutive law."""
        return {'E': self.Ec}

    def __parabolarectangle__(self) -> dict:
        """Returns kwargs for creating a parabola-rectangle constitutive law.

        Uses Hognestad peak strain (0.002) and ACI ultimate strain (0.003).
        """
        return {
            'fc': self.fcd(),
            'eps_0': 0.002,
            'eps_u': self.eps_cu,
            'n': 2,
        }

    def __bilinearcompression__(self) -> dict:
        """Returns kwargs for creating a bilinear compression law."""
        return {
            'fc': self.fcd(),
            'eps_c': 0.002,
            'eps_cu': self.eps_cu,
        }
```

- [ ] **Step 4: Register in the factory**

In `structuralcodes/materials/concrete/__init__.py`, add the import:

```python
from ._concreteACI318_25 import ConcreteACI318_25
```

Add to `__all__`:
```python
'ConcreteACI318_25',
```

Add to `CONCRETES`:
```python
CONCRETES: t.Dict[str, Concrete] = {
    'ACI 318-25': ConcreteACI318_25,
    'fib Model Code 2010': ConcreteMC2010,
    'EUROCODE 2 1992-1-1:2004': ConcreteEC2_2004,
    'EUROCODE 2 1992-1-1:2023': ConcreteEC2_2023,
}
```

Note: The key `'ACI 318-25'` must match `__title__` in the code module's `__init__.py`.

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_aci318_25/test_concrete_aci318_25.py -v`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add structuralcodes/materials/concrete/_concreteACI318_25.py structuralcodes/materials/concrete/__init__.py tests/test_aci318_25/test_concrete_aci318_25.py
git commit -m "feat(aci318_25): add ConcreteACI318_25 material class"
```

---

## Task 9: ReinforcementACI318_25 Material Class

**Files:**
- Create: `structuralcodes/materials/reinforcement/_reinforcementACI318_25.py`
- Create: `tests/test_aci318_25/test_reinforcement_aci318_25.py`
- Modify: `structuralcodes/materials/reinforcement/__init__.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_aci318_25/test_reinforcement_aci318_25.py
"""Tests for the ReinforcementACI318_25 material class."""

import math

import pytest

import structuralcodes
from structuralcodes.materials.reinforcement._reinforcementACI318_25 import (
    ReinforcementACI318_25,
)


@pytest.fixture(autouse=True)
def _reset_design_code():
    yield
    structuralcodes.set_design_code(None)


def _make_gr60(**kwargs):
    return ReinforcementACI318_25(
        fyk=420, Es=200000, ftk=550, epsuk=0.05, **kwargs,
    )


class TestConstruction:
    def test_basic(self):
        r = _make_gr60()
        assert r.fyk == 420

    def test_default_name(self):
        r = _make_gr60()
        assert r.name == 'Reinforcement420'

    def test_from_grade(self):
        r = ReinforcementACI318_25.from_grade('60')
        assert r.fyk == 420
        assert r.ftk == 550


class TestProperties:
    def test_gamma_s(self):
        r = _make_gr60()
        assert r.gamma_s == 1.0

    def test_fyd(self):
        r = _make_gr60()
        assert math.isclose(r.fyd(), 420)

    def test_ftd(self):
        r = _make_gr60()
        assert math.isclose(r.ftd(), 550)

    def test_epsud(self):
        r = _make_gr60()
        assert r.epsud() == 0.05


class TestConstitutiveLaws:
    def test_elastic(self):
        r = _make_gr60(constitutive_law='elastic')
        assert r.constitutive_law is not None

    def test_elasticplastic(self):
        r = _make_gr60(constitutive_law='elasticplastic')
        assert r.constitutive_law is not None

    def test_elasticperfectlyplastic(self):
        r = _make_gr60(constitutive_law='elasticperfectlyplastic')
        assert r.constitutive_law is not None


class TestFactory:
    def test_create_via_factory(self):
        from structuralcodes.materials.reinforcement import create_reinforcement
        r = create_reinforcement(
            fyk=420, Es=200000, ftk=550, epsuk=0.05, design_code='aci318_25',
        )
        assert isinstance(r, ReinforcementACI318_25)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_aci318_25/test_reinforcement_aci318_25.py -v`

Expected: ERRORS — module does not exist.

- [ ] **Step 3: Implement ReinforcementACI318_25**

```python
# structuralcodes/materials/reinforcement/_reinforcementACI318_25.py
"""Reinforcement material class for ACI 318-25."""

from __future__ import annotations

import typing as t

from structuralcodes.codes import aci318_25

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._reinforcement import Reinforcement


class ReinforcementACI318_25(Reinforcement):
    """ACI 318-25 reinforcement material.

    Strengths are unreduced (gamma_s=1.0). ACI applies strength reduction
    factors at member capacity level, not material level.

    Supports ASTM A615 grades: 40, 60, 80, 100.
    """

    def __init__(
        self,
        fyk: float,
        Es: float = 200000.0,
        ftk: float = 550.0,
        epsuk: float = 0.05,
        gamma_s: t.Optional[float] = None,
        name: t.Optional[str] = None,
        density: float = 7850,
        constitutive_law: t.Optional[
            t.Union[
                t.Literal['elastic', 'elasticplastic', 'elasticperfectlyplastic'],
                ConstitutiveLaw,
            ]
        ] = 'elasticperfectlyplastic',
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[bool] = None,
        **kwargs,
    ) -> None:
        """Initialize ACI 318-25 reinforcement.

        Args:
            fyk: Specified yield strength fy in MPa.
            Es: Modulus of elasticity in MPa (default 200000).
            ftk: Ultimate tensile strength fu in MPa (default 550).
            epsuk: Ultimate strain (default 0.05).
            gamma_s: Partial factor (default 1.0 for ACI, no material reduction).
            name: Descriptive name.
            density: Density in kg/m3 (default 7850).
            constitutive_law: ConstitutiveLaw or string name.
            initial_strain: Initial strain of the material.
            initial_stress: Initial stress of the material.
            strain_compatibility: Whether material deforms with geometry.
        """
        del kwargs
        if name is None:
            name = f'Reinforcement{round(fyk):d}'

        super().__init__(
            fyk=fyk,
            Es=Es,
            name=name,
            density=density,
            ftk=ftk,
            epsuk=epsuk,
            gamma_s=gamma_s,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            strain_compatibility=strain_compatibility,
        )
        self._constitutive_law = (
            constitutive_law
            if isinstance(constitutive_law, ConstitutiveLaw)
            else create_constitutive_law(
                constitutive_law_name=constitutive_law, material=self
            )
        )
        if 'steel' not in self._constitutive_law.__materials__:
            raise ValueError(
                'The provided constitutive law is not valid for reinforcement.'
            )
        self._apply_initial_strain()

    @classmethod
    def from_grade(
        cls,
        grade: str = '60',
        epsuk: float = 0.05,
        **kwargs,
    ) -> 'ReinforcementACI318_25':
        """Create from ASTM A615 grade.

        Args:
            grade: ASTM grade string ('40', '60', '80', '100').
            epsuk: Ultimate strain (default 0.05).

        Returns:
            ReinforcementACI318_25 instance.
        """
        props = aci318_25.reinforcement_grade_props(grade)
        return cls(
            fyk=props['fy'],
            ftk=props['fu'],
            epsuk=epsuk,
            **kwargs,
        )

    def fyd(self) -> float:
        """Design yield strength.

        ACI 318 does not reduce material strength. Returns fy / gamma_s,
        which with default gamma_s=1.0 gives the unreduced yield strength.
        """
        return self.fyk / self.gamma_s

    @property
    def gamma_s(self) -> float:
        """Partial factor for reinforcement.

        Default is 1.0 for ACI 318 (no material partial factor).
        """
        return self._gamma_s or 1.0

    def ftd(self) -> float:
        """Design ultimate strength."""
        return self.ftk / self.gamma_s

    def epsud(self) -> float:
        """Design ultimate strain."""
        return self.epsuk

    def __elastic__(self) -> dict:
        """Returns kwargs for an elastic constitutive law."""
        return {'E': self.Es}

    def __elasticperfectlyplastic__(self) -> dict:
        """Returns kwargs for ElasticPlastic law with no hardening."""
        return {
            'E': self.Es,
            'fy': self.fyd(),
            'eps_su': self.epsud(),
        }

    def __elasticplastic__(self) -> dict:
        """Returns kwargs for ElasticPlastic law with hardening."""
        Eh = (self.ftd() - self.fyd()) / (self.epsud() - self.epsyd)
        return {
            'E': self.Es,
            'fy': self.fyd(),
            'Eh': Eh,
            'eps_su': self.epsud(),
        }
```

- [ ] **Step 4: Register in the factory**

In `structuralcodes/materials/reinforcement/__init__.py`, add:

```python
from ._reinforcementACI318_25 import ReinforcementACI318_25
```

Add to `__all__`:
```python
'ReinforcementACI318_25',
```

Add to `REINFORCEMENTS`:
```python
REINFORCEMENTS: t.Dict[str, Reinforcement] = {
    'ACI 318-25': ReinforcementACI318_25,
    'fib Model Code 2010': ReinforcementMC2010,
    'EUROCODE 2 1992-1-1:2004': ReinforcementEC2_2004,
    'EUROCODE 2 1992-1-1:2023': ReinforcementEC2_2023,
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_aci318_25/test_reinforcement_aci318_25.py -v`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add structuralcodes/materials/reinforcement/_reinforcementACI318_25.py structuralcodes/materials/reinforcement/__init__.py tests/test_aci318_25/test_reinforcement_aci318_25.py
git commit -m "feat(aci318_25): add ReinforcementACI318_25 material class"
```

---

## Task 10: Whitney Block Constitutive Law

**Files:**
- Create: `structuralcodes/materials/constitutive_laws/_whitneyblock.py`
- Create: `tests/test_aci318_25/test_whitneyblock.py`
- Modify: `structuralcodes/materials/constitutive_laws/__init__.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_aci318_25/test_whitneyblock.py
"""Tests for the WhitneyBlock constitutive law."""

import math

import numpy as np
import pytest

from structuralcodes.materials.constitutive_laws._whitneyblock import WhitneyBlock


@pytest.fixture
def wb():
    """Whitney block for 4000 psi concrete: fc=0.85*27.58=23.44, beta1=0.85."""
    return WhitneyBlock(fc=23.44, beta1=0.85, eps_cu=0.003)


class TestGetStress:
    def test_in_active_zone(self, wb):
        """Strain in the active zone should return -fc."""
        # Active zone: eps between -0.003 and -0.003*(1-0.85) = -0.00045
        eps = -0.002  # well within active zone
        stress = wb.get_stress(eps)
        assert math.isclose(stress, -23.44, rel_tol=1e-6)

    def test_in_zero_zone(self, wb):
        """Strain below the active zone should return 0."""
        eps = -0.0002  # between 0 and -0.00045
        stress = wb.get_stress(eps)
        assert stress == 0.0

    def test_positive_strain(self, wb):
        """Positive (tensile) strain returns 0."""
        assert wb.get_stress(0.001) == 0.0

    def test_beyond_ultimate(self, wb):
        """Strain beyond eps_cu returns 0."""
        assert wb.get_stress(-0.004) == 0.0

    def test_array_input(self, wb):
        eps = np.array([-0.004, -0.002, -0.0002, 0.001])
        sig = wb.get_stress(eps)
        expected = np.array([0.0, -23.44, 0.0, 0.0])
        np.testing.assert_allclose(sig, expected, atol=1e-6)


class TestGetUltimateStrain:
    def test_values(self, wb):
        eps_min, eps_max = wb.get_ultimate_strain()
        assert math.isclose(eps_min, -0.003)
        assert eps_max == 0.0


class TestGetTangent:
    def test_in_active_zone(self, wb):
        """Tangent is 0 everywhere (piecewise constant)."""
        assert wb.get_tangent(-0.002) == 0.0

    def test_at_zero(self, wb):
        assert wb.get_tangent(0.0) == 0.0


class TestMarin:
    def test_returns_strains_and_coefficients(self, wb):
        """Marin integration should return valid strain limits and coefficients."""
        strain = (-0.0015, -0.001)  # linear strain profile
        strains, coeff = wb.__marin__(strain)
        assert strains is not None or coeff is not None
        # Should have at least one region
        assert len(coeff) >= 1
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_aci318_25/test_whitneyblock.py -v`

Expected: ERRORS — module does not exist.

- [ ] **Step 3: Implement WhitneyBlock**

```python
# structuralcodes/materials/constitutive_laws/_whitneyblock.py
"""Whitney equivalent rectangular stress block constitutive law."""

from __future__ import annotations

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


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
            Stored internally as a negative value (compression).
        beta1: Depth factor mapping neutral axis depth c to block
            depth a = beta1*c.
        eps_cu: Ultimate concrete strain (default 0.003).
    """

    __materials__: t.Tuple[str, ...] = ('concrete',)

    def __init__(
        self,
        fc: float,
        beta1: float,
        eps_cu: float = 0.003,
        name: t.Optional[str] = None,
    ) -> None:
        name = name if name is not None else 'WhitneyBlock'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._beta1 = beta1
        self._eps_cu = -abs(eps_cu)
        self._eps_transition = self._eps_cu * (1.0 - beta1)

    def get_stress(
        self, eps: t.Union[float, ArrayLike],
    ) -> t.Union[float, ArrayLike]:
        """Return stress for given strain.

        Returns -fc (compression) in the active zone, 0 elsewhere.
        """
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        eps = self.preprocess_strains_with_limits(eps=eps)

        if np.isscalar(eps):
            if self._eps_cu <= eps <= self._eps_transition:
                return self._fc
            return 0.0

        sig = np.zeros_like(eps, dtype=float)
        active = (eps >= self._eps_cu) & (eps <= self._eps_transition)
        sig[active] = self._fc
        return sig

    def get_tangent(
        self, eps: t.Union[float, ArrayLike],
    ) -> t.Union[float, ArrayLike]:
        """Return tangent modulus. Always 0 (piecewise constant)."""
        if np.isscalar(eps):
            return 0.0
        return np.zeros_like(np.atleast_1d(eps), dtype=float)

    def get_ultimate_strain(
        self, yielding: bool = False,
    ) -> t.Tuple[float, float]:
        """Return ultimate strain (negative, positive)."""
        return (self._eps_cu, 0.0)

    def __marin__(
        self, strain: t.Tuple[float, float],
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
        """Returns strain limits and coefficients for Marin integration.

        The Whitney block has two regions:
        - Zero stress zone: from eps_transition to 0
        - Constant stress zone: from eps_cu to eps_transition

        Args:
            strain: Tuple (eps_0, eps_1) defining linear strain profile.

        Returns:
            (strains, coeff) for Marin integration.
        """
        strains = []
        coeff = []

        if strain[1] == 0:
            # Uniform strain
            eps_0 = self.preprocess_strains_with_limits(strain[0])
            if self._eps_cu <= eps_0 <= self._eps_transition:
                strains = None
                coeff.append((self._fc,))
            else:
                strains = None
                coeff.append((0.0,))
        else:
            # Constant stress zone
            strains.append((self._eps_cu, self._eps_transition))
            coeff.append((self._fc,))
            # Zero stress zone
            strains.append((self._eps_transition, 0))
            coeff.append((0.0,))

        return strains, coeff

    def __marin_tangent__(
        self, strain: t.Tuple[float, float],
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
        """Returns strain limits and coefficients for Marin tangent integration.

        Tangent is always 0 for piecewise constant law.
        """
        strains = []
        coeff = []

        if strain[1] == 0:
            strains = None
            coeff.append((0.0,))
        else:
            strains.append((self._eps_cu, 0))
            coeff.append((0.0,))

        return strains, coeff
```

- [ ] **Step 4: Register in the constitutive laws factory**

In `structuralcodes/materials/constitutive_laws/__init__.py`, add the import:

```python
from ._whitneyblock import WhitneyBlock
```

Add `'WhitneyBlock'` to `__all__`.

Add to `CONSTITUTIVE_LAWS`:
```python
'whitneyblock': WhitneyBlock,
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_aci318_25/test_whitneyblock.py -v`

Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add structuralcodes/materials/constitutive_laws/_whitneyblock.py structuralcodes/materials/constitutive_laws/__init__.py tests/test_aci318_25/test_whitneyblock.py
git commit -m "feat: add WhitneyBlock constitutive law for equivalent rectangular stress block"
```

---

## Task 11: End-to-End One-Way Slab Example

**Files:**
- Create: `tests/test_aci318_25/test_one_way_slab_example.py`

- [ ] **Step 1: Write the end-to-end test**

```python
# tests/test_aci318_25/test_one_way_slab_example.py
"""End-to-end validation: one-way slab design with both paths.

Problem: 4000 psi concrete, Gr 60 rebar, 20 ft span, one end continuous.
Design per 12 in. strip.
"""

import math

import pytest
from shapely.geometry import Polygon

import structuralcodes
from structuralcodes.codes import aci318_25
from structuralcodes.geometry import (
    CompoundGeometry,
    PointGeometry,
    SurfaceGeometry,
)
from structuralcodes.materials.concrete import ConcreteACI318_25
from structuralcodes.materials.reinforcement import ReinforcementACI318_25
from structuralcodes.sections import BeamSection

# Problem parameters (SI)
FC = 27.58     # MPa (4000 psi)
FY = 420.0     # MPa (60 ksi)
SPAN = 6096.0  # mm (20 ft)
B = 305.0      # mm (12 in. strip)


@pytest.fixture(autouse=True)
def _reset_design_code():
    yield
    structuralcodes.set_design_code(None)


class TestPathA_ClosedForm:
    """Path A: closed-form ACI equations."""

    def test_min_thickness(self):
        h = aci318_25.min_thickness(SPAN, 'one_end_continuous')
        assert math.isclose(h, SPAN / 24, rel_tol=1e-6)
        # Round up to practical thickness
        assert h > 200  # should be ~254 mm (10 in.)

    def test_flexure_design(self):
        h = 254.0  # 10 in.
        d = h - 19 - 16 / 2  # 3/4" cover + half #5 bar = 227 mm

        # Assume Mu = 40e6 N-mm for this example
        Mu = 40e6
        phi = 0.9
        As = aci318_25.As_required(Mu, phi, FY, FC, B, d)
        As_min = aci318_25.As_min_slab(FY, B, h)
        As_design = max(As, As_min)

        # Check tension-controlled
        a = aci318_25.stress_block_depth_sr(As_design, FY, FC, B)
        c = aci318_25.neutral_axis_depth(a, aci318_25.beta1(FC))
        eps_t = aci318_25.eps_t_from_c(c, d)
        assert aci318_25.As_max_check(eps_t, FY)
        assert aci318_25.phi_flexure(eps_t, FY) == 0.9

        # Verify Mn
        Mn = aci318_25.Mn_singly_reinforced(As_design, FY, FC, B, d)
        assert phi * Mn >= Mu

    def test_shear_check(self):
        h = 254.0
        d = 227.0
        rho_w = 0.009  # typical

        Vc = aci318_25.Vc_detailed(FC, B, d, rho_w)
        phi_Vc = aci318_25.phi_shear() * Vc

        # For a typical slab, Vu should be less than phi*Vc
        # Use a reasonable Vu (e.g., 30 kN)
        Vu = 30000.0  # N
        assert not aci318_25.shear_reinforcement_required(Vu, phi_Vc)


class TestPathB_SectionIntegrator:
    """Path B: section integrator with ACI materials."""

    def test_section_analysis(self):
        concrete = ConcreteACI318_25(fck=FC, constitutive_law='parabolarectangle')
        steel = ReinforcementACI318_25(
            fyk=FY, Es=200000, ftk=550, epsuk=0.05,
            constitutive_law='elasticperfectlyplastic',
        )

        poly = Polygon([(0, 0), (B, 0), (B, 254), (0, 254)])
        surf = SurfaceGeometry(poly, concrete)
        bar1 = PointGeometry(point=(100, 27), diameter=16, material=steel)
        bar2 = PointGeometry(point=(205, 27), diameter=16, material=steel)
        section_geo = CompoundGeometry([surf], [bar1, bar2])

        section = BeamSection(section_geo, integrator='marin')
        props = section.gross_properties
        assert props.area > 0

    def test_factory_round_trip(self):
        """Verify material factory works with aci318_25 code."""
        from structuralcodes.materials.concrete import create_concrete
        from structuralcodes.materials.reinforcement import create_reinforcement

        structuralcodes.set_design_code('aci318_25')
        c = create_concrete(fck=FC)
        assert isinstance(c, ConcreteACI318_25)

        r = create_reinforcement(fyk=FY, Es=200000, ftk=550, epsuk=0.05)
        assert isinstance(r, ReinforcementACI318_25)


class TestCrossCheck:
    """Cross-check between Path A and Path B."""

    def test_mn_agreement(self):
        """Closed-form Mn and integrator Mn should agree within 5%."""
        As = 400.0  # mm2 (2 x #5 bars)
        d = 227.0
        h = 254.0

        # Path A: closed-form
        Mn_closed = aci318_25.Mn_singly_reinforced(As, FY, FC, B, d)

        # Path B: integrator
        concrete = ConcreteACI318_25(fck=FC, constitutive_law='parabolarectangle')
        steel = ReinforcementACI318_25(
            fyk=FY, Es=200000, ftk=550, epsuk=0.05,
            constitutive_law='elasticperfectlyplastic',
        )

        poly = Polygon([(0, 0), (B, 0), (B, h), (0, h)])
        surf = SurfaceGeometry(poly, concrete)
        # Place bars to match As=400 mm2 (2 bars of ~200 mm2 each, dia ~16)
        bar1 = PointGeometry(point=(100, 27), diameter=16, material=steel)
        bar2 = PointGeometry(point=(205, 27), diameter=16, material=steel)
        section_geo = CompoundGeometry([surf], [bar1, bar2])

        section = BeamSection(section_geo, integrator='marin')
        calc = section.section_calculator
        strain = calc.find_equilibrium_fixed_pivot(
            geom=section.geometry, n=0, yielding=True,
        )
        N, My, Mz, data = calc.integrator.integrate_strain_response_on_geometry(
            geometry=section.geometry, strain=strain,
        )
        Mn_integrator = abs(My)

        # Allow 5% tolerance (Whitney vs parabola-rectangle)
        assert math.isclose(Mn_closed, Mn_integrator, rel_tol=0.05)
```

- [ ] **Step 2: Run the end-to-end tests**

Run: `pytest tests/test_aci318_25/test_one_way_slab_example.py -v`

Expected: All tests PASS.

- [ ] **Step 3: Run the full test suite to check for regressions**

Run: `pytest --tb=short -q`

Expected: All existing tests still pass. New ACI tests pass. Zero regressions.

- [ ] **Step 4: Commit**

```bash
git add tests/test_aci318_25/test_one_way_slab_example.py
git commit -m "test(aci318_25): add end-to-end one-way slab validation (both paths)"
```

---

## Task 12: Final Regression Check and Cleanup

- [ ] **Step 1: Run ruff formatting**

Run: `ruff format structuralcodes/codes/aci318_25/ structuralcodes/materials/concrete/_concreteACI318_25.py structuralcodes/materials/reinforcement/_reinforcementACI318_25.py structuralcodes/materials/constitutive_laws/_whitneyblock.py tests/test_aci318_25/`

- [ ] **Step 2: Run ruff linting**

Run: `ruff check structuralcodes/codes/aci318_25/ structuralcodes/materials/concrete/_concreteACI318_25.py structuralcodes/materials/reinforcement/_reinforcementACI318_25.py structuralcodes/materials/constitutive_laws/_whitneyblock.py tests/test_aci318_25/`

Fix any issues found.

- [ ] **Step 3: Run full test suite**

Run: `pytest --tb=short -q`

Expected: All tests pass, zero regressions.

- [ ] **Step 4: Commit any formatting/lint fixes**

```bash
git add -u
git commit -m "style(aci318_25): fix formatting and lint issues"
```
