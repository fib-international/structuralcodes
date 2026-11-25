"""Tests for the functions in corrosion."""

from math import isclose

import numpy as np
import pytest
from scipy.stats import norm

from structuralcodes.codes.mc2020._corrosion import (
    kc,
    pcorr_fractile,
    pcorr_rep,
)


@pytest.mark.parametrize(
    'corrosion_type, exposure_class, expected',
    [
        ('carbonation', 'sheltered', (2, 3)),
        ('carbonation', 'unsheltered', (5, 7)),
        ('chloride', 'wet', (4, 6)),
        ('chloride', 'cyclic_dry_wet', (30, 40)),
        ('chloride', 'airborn_seawater', (30, 40)),
        ('chloride', 'submerged', (4, 7)),
        ('chloride', 'tidal_zone', (50, 100)),
        ('Chloride', 'Tidal_zone', (50, 100)),
    ],
)
def test_pcorr_dict(corrosion_type, exposure_class, expected):
    """Test pcorr_rep with valid input."""
    pcorr_dic = pcorr_rep(corrosion_type, exposure_class)
    assert pcorr_dic['mean'] == expected[0]
    assert pcorr_dic['std'] == expected[1]


@pytest.mark.parametrize(
    'corrosion_type, exposure_class',
    [
        ('chloride', 'sheltered'),
        ('carbonatio', 'unsheltered'),
        ('chloride', 'wett'),
        ('chloride', 'sheltered'),
        ('chloride', 'sunmerged'),
    ],
)
def test_pcorr_dict_invalid(corrosion_type, exposure_class):
    """Test pcorr_rep with valid input."""
    with pytest.raises(ValueError):
        pcorr_rep(corrosion_type, exposure_class)


@pytest.mark.parametrize(
    'fractile',
    [
        0.5,
        0.05,
        0.95,
        0.16,
        0.84,
    ],
)
@pytest.mark.parametrize(
    'corrosion_type, exposure_class',
    [
        ('carbonation', 'sheltered'),
        ('carbonation', 'unsheltered'),
        ('chloride', 'wet'),
        ('chloride', 'cyclic_dry_wet'),
        ('chloride', 'airborn_seawater'),
        ('chloride', 'submerged'),
        ('chloride', 'tidal_zone'),
    ],
)
def test_pcorr_fractile(corrosion_type, exposure_class, fractile):
    """Test pcorr_fractile with valid input."""
    # Evaluate expected value
    pcorr_dic = pcorr_rep(
        corrosion_type=corrosion_type, exposure_class=exposure_class
    )
    scale = (
        pcorr_dic['mean'] ** 2
        / (pcorr_dic['mean'] ** 2 + pcorr_dic['std'] ** 2) ** 0.5
    )
    median = scale
    mu = np.log(median)  # median = scale = exp(mu)
    s = (
        np.log(1 + pcorr_dic['std'] ** 2 / pcorr_dic['mean'] ** 2) ** 0.5
    )  # sigma = s
    z = norm.ppf(fractile)
    expected = np.exp(mu + z * s)
    # Act
    assert isclose(
        pcorr_fractile(
            corrosion_type=corrosion_type,
            exposure_class=exposure_class,
            fractile=fractile,
        ),
        expected,
    )


@pytest.mark.parametrize(
    'corrosion_type, exposure_class',
    [
        ('chloride', 'sheltered'),
        ('carbonatio', 'unsheltered'),
        ('chloride', 'wett'),
        ('chloride', 'sheltered'),
        ('chloride', 'sunmerged'),
    ],
)
def test_pcorr_fractile_invalid(corrosion_type, exposure_class):
    """Test pcorr_fractile with invalid  corrosion type and/or exposure."""
    with pytest.raises(ValueError):
        (
            pcorr_fractile(
                corrosion_type=corrosion_type,
                exposure_class=exposure_class,
                fractile=0.5,
            ),
        )


@pytest.mark.parametrize(
    'fractile',
    [
        0.0,
        -0.1,
        50,
        1.0,
        1.2,
        95.0,
    ],
)
@pytest.mark.parametrize(
    'corrosion_type, exposure_class',
    [
        ('carbonation', 'sheltered'),
        ('carbonation', 'unsheltered'),
        ('chloride', 'wet'),
        ('chloride', 'cyclic_dry_wet'),
        ('chloride', 'airborn_seawater'),
        ('chloride', 'submerged'),
        ('chloride', 'tidal_zone'),
    ],
)
def test_pcorr_fractile_invalid_fractile(
    corrosion_type, exposure_class, fractile
):
    """Test pcorr_fractile with invalid fractile."""
    with pytest.raises(ValueError):
        (
            pcorr_fractile(
                corrosion_type=corrosion_type,
                exposure_class=exposure_class,
                fractile=fractile,
            ),
        )


@pytest.mark.parametrize(
    'fc, expected',
    [
        (20, 0.75),
        (25, 0.75),
        (30, 0.75),
        (35, 0.712435688694747),
        (40, 0.681420222312052),
        (45, 0.655185348552224),
        (50, 0.632574498976312),
    ],
)
def test_kc(fc, expected):
    """Test kc with valid input."""
    assert isclose(kc(fc), expected)


@pytest.mark.parametrize(
    'fc',
    [
        (0),
        (-20),
    ],
)
def test_kc_invalid(fc):
    """Test kc with invalid input."""
    with pytest.raises(ValueError):
        kc(fc)
