"""Tests of functions from Section 6 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section6_durability


@pytest.mark.parametrize(
    'env_level, env_name, expected_classes',
    [
        (
            None,
            None,
            [
                'X0',
                'XC1',
                'XC2',
                'XC3',
                'XC4',
                'XD1',
                'XD2',
                'XD3',
                'XS1',
                'XS2',
                'XS3',
                'XF1',
                'XF2',
                'XF3',
                'XF4',
                'XA1',
                'XA2',
                'XA3',
                'XM1',
                'XM2',
                'XM3',
            ],
        ),
        (1, None, ['X0']),
        (2, None, ['XC1', 'XC2', 'XC3', 'XC4']),
        (3, None, ['XD1', 'XD2', 'XD3']),
        (4, None, ['XS1', 'XS2', 'XS3']),
        (5, None, ['XF1', 'XF2', 'XF3', 'XF4']),
        (6, None, ['XA1', 'XA2', 'XA3']),
        (7, None, ['XM1', 'XM2', 'XM3']),
        (None, 'none', ['X0']),
        (None, 'carbonation', ['XC1', 'XC2', 'XC3', 'XC4']),
        (None, 'chlorides', ['XD1', 'XD2', 'XD3']),
        (None, 'sea', ['XS1', 'XS2', 'XS3']),
        (None, 'freeze', ['XF1', 'XF2', 'XF3', 'XF4']),
        (None, 'chemical', ['XA1', 'XA2', 'XA3']),
        (None, 'abrasion', ['XM1', 'XM2', 'XM3']),
    ],
)
def test_get_exposure_classes(env_level, env_name, expected_classes):
    """Test if the function returns the correct list of exposure classes."""
    assert (
        _section6_durability.get_exposure_classes(env_level, env_name)
        == expected_classes
    )


@pytest.mark.parametrize(
    'exp_class, expected_main_description',
    [
        ('X0', 'No risk of corrosion or attack.'),
        ('XC1', 'Corrosion of embedded metal induced by carbonation.'),
        (
            'XD3',
            'Corrosion of embedded metal induced by chlorides, '
            + 'excluding sea water.',
        ),
    ],
)
def test_get_exposure_class_description_valid(
    exp_class, expected_main_description
):
    """Test valid exposure classes."""
    description = _section6_durability.get_exposure_class_description(
        exp_class
    )
    assert description['main_description'] == expected_main_description


def test_get_exporuse_class_description_invalid():
    """Test invalid exposure class."""
    with pytest.raises(ValueError):
        _section6_durability.get_exposure_class_description('INVALID_CLASS')

    # Test another invalid input
    with pytest.raises(ValueError):
        _section6_durability.get_exposure_class_description('UNKNOWN')


def test_get_exporuse_class_description_case_insensitive():
    """Test if the function handles case-insensitive class codes."""
    assert (
        'Corrosion of embedded metal induced by carbonation.'
        in _section6_durability.get_exposure_class_description('xc1')[
            'main_description'
        ]
    )
    assert (
        'Corrosion of embedded metal induced by chlorides'
        in _section6_durability.get_exposure_class_description('xd3')[
            'main_description'
        ]
    )


def test_cnom():
    """Test the nominal cover calculation."""
    assert _section6_durability.c_nom(25, 5) == 30
    assert _section6_durability.c_nom(20, 10) == 30
    assert _section6_durability.c_nom(0, 0) == 0

    with pytest.raises(ValueError):
        _section6_durability.c_nom(-1, 5)

    with pytest.raises(ValueError):
        _section6_durability.c_nom(25, -5)


def test_cmin():
    """Test the minimum cover calculation."""
    assert _section6_durability.c_min(25, 5, 20) == 30
    assert _section6_durability.c_min(20, 10, 20) == 30
    assert _section6_durability.c_min(10, 0, 15) == 15
    assert _section6_durability.c_min(5, 0, 5, 5) == 15

    with pytest.raises(ValueError):
        _section6_durability.c_min(-1, 5, 20)

    with pytest.raises(ValueError):
        _section6_durability.c_min(25, 5, -20)

    with pytest.raises(ValueError):
        _section6_durability.c_min(25, 5, 20, -5)


@pytest.mark.parametrize(
    'exposure_class, design_service_life, exposure_resistance_class, expected',
    [
        ('XC1', 50, 'XRC 0.5', 10),
        ('XC3', 100, 'XRC 3', 30),
        ('XC3', 50, 'XRC 2', 15),
        ('XC4', 100, 'XRC 4', 40),
        ('XC2', 75, 'XRC 7', 32.5),
    ],
)
def test_minimum_cover_carbonation(
    exposure_class, design_service_life, exposure_resistance_class, expected
):
    """Test minimum concrete cover for carbonation."""
    assert (
        _section6_durability.c_min_dur_carb(
            exposure_class, design_service_life, exposure_resistance_class
        )
        == expected
    )


@pytest.mark.parametrize(
    'exposure_class, design_service_life, exposure_resistance_class, expected',
    [
        ('XS1', 50, 'XRDS 0.5', 20),
        ('XS3', 100, 'XRDS 3', 65),
        ('XD1', 50, 'XRDS 2', 25),
        ('XD3', 100, 'XRDS 1', 45),
        ('XD3', 75, 'XRDS 1', 40),
        ('XD2', 50, 'XRDS 10', 80),
    ],
)
def test_minimum_cover_chlorides(
    exposure_class, design_service_life, exposure_resistance_class, expected
):
    """Test minimum concrete cover for chlorides."""
    assert (
        _section6_durability.c_min_dur_chlo(
            exposure_class, design_service_life, exposure_resistance_class
        )
        == expected
    )


@pytest.mark.parametrize(
    'exposure_class, design_service_life, exposure_resistance_class',
    [
        ('XS3', 100, 'XRDS 10'),
        ('XD3', 100, 'XRDS 8'),
    ],
)
def test_minimum_cover_chlorides_with_invalid_values(
    exposure_class, design_service_life, exposure_resistance_class
):
    """Test minimum concrete cover for chlorides."""
    with pytest.raises(ValueError):
        _section6_durability.c_min_dur_chlo(
            exposure_class, design_service_life, exposure_resistance_class
        )


def test_delta_c_min_30():
    """Test reduction of minimum cover for design life of 30 years."""
    assert _section6_durability.delta_c_min_30() == -5.0


def test_delta_c_min_exc():
    """Test reduction of minimum cover for superior
    compaction or improved curing.
    """
    assert _section6_durability.delta_c_min_exc() == -5.0


def test_delta_c_min_p():
    """Test addition of minimum cover for prestressing tendons."""
    assert _section6_durability.delta_c_min_p() == 10.0


def test_delta_dur_red_1():
    """Test addition of minimum cover for special
    reinforcing steel measures.
    """
    assert _section6_durability.delta_dur_red_1() == -10.0


def test_delta_dur_red_2():
    """Test reduction of minimum cover for
    special reinforcing steel measures.
    """
    assert _section6_durability.delta_dur_red_2() == 0.0


@pytest.mark.parametrize(
    'xm_class, expected',
    [
        ('XM1', 5.0),
        ('XM2', 10.0),
        ('XM3', 15.0),
    ],
)
def test_delta_dur_abr(xm_class, expected):
    """Test addition of minimum cover for abrasion based on XM class."""
    assert (
        _section6_durability.delta_dur_abr(xm_class) == expected
    ), f'Expected addition of {expected} mm for XM class {xm_class}.'
