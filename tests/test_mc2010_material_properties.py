"""Tests for the _concrete_material_properties module"""
import math
import numpy as np

import pytest

from structuralcodes.codes.mc2010 import _concrete_material_properties


@pytest.mark.parametrize(
    'test_input, expected',
    [(12, 20), (35, 43), (55, 63), (90, 98), (120, 128)],
)
def test_fcm(test_input, expected):
    """Test the fcm function."""
    assert math.isclose(
        _concrete_material_properties.fcm(test_input), expected
    )


@pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 1.6),
        (16, 1.9),
        (20, 2.2),
        (25, 2.6),
        (30, 2.9),
        (35, 3.2),
        (40, 3.5),
        (45, 3.8),
        (50, 4.1),
        (55, 4.2),
        (60, 4.4),
        (70, 4.6),
        (80, 4.8),
        (90, 5.0),
        (100, 5.2),
        (110, 5.4),
        (120, 5.6),
    ],
)
def test_fctm(test_input, expected):
    """Test the fctm function."""
    assert math.isclose(
        _concrete_material_properties.fctm(test_input), expected, rel_tol=0.02
    )


@pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 1.1),
        (16, 1.3),
        (20, 1.5),
        (25, 1.8),
        (30, 2.0),
        (35, 2.2),
        (40, 2.5),
        (45, 2.7),
        (50, 2.9),
        (55, 3.0),
        (60, 3.1),
        (70, 3.2),
        (80, 3.4),
        (90, 3.5),
        (100, 3.7),
        (110, 3.8),
        (120, 3.9),
    ],
)
def test_fctkmin(test_input, expected):
    """Test the fctkmin function."""
    assert math.isclose(
        _concrete_material_properties.fctkmin(
            _concrete_material_properties.fctm(test_input)
        ),
        expected,
        rel_tol=0.031,
    )


@pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 2.0),
        (16, 2.5),
        (20, 2.9),
        (25, 3.3),
        (30, 3.8),
        (35, 4.2),
        (40, 4.6),
        (45, 4.8),
        (50, 5.3),
        (55, 5.5),
        (60, 5.7),
        (70, 6.0),
        (80, 6.3),
        (90, 6.6),
        (100, 6.8),
        (110, 7.0),
        (120, 7.2),
    ],
)
def test_fctkmax(test_input, expected):
    """Test the fctkmax function."""
    assert math.isclose(
        _concrete_material_properties.fctkmax(
            _concrete_material_properties.fctm(test_input)
        ),
        expected,
        rel_tol=0.028,
    )


@pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 125.172),
        (35, 143.664),
        (55, 153.888),
        (90, 166.626),
        (120, 174.832),
    ],
)
def test_Gf(test_input, expected):
    """Test the Gf function."""
    assert math.isclose(
        _concrete_material_properties.Gf(test_input),
        expected,
        rel_tol=1e-5,
    )


@pytest.mark.parametrize(
    'agg_type, expected',
    [
        ('basalt', None),
        ('quartzite', None),
        ('Limestone', None),
        ('SANDSTONE', None),
        ('test', None) # Raises ValueError
    ],
)
def test_check_agg_type(agg_type, expected):
    """Test check_agg_type function."""
    try:
        assert _concrete_material_properties._check_agg_types(agg_type) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            _concrete_material_properties._check_agg_types(agg_type) == expected
            assert str(exc_info.value).startswith(
                    "The specified cement type")


@pytest.mark.parametrize(
    'cem_class, expected',
    [
        ("32.5 N", None),
        ("32.5 R", None),
        ("42.5 N", None),
        ("42.5 R", None),
        ("52.5 N", None),
        ("52.5 R", None),
        ("", None), # raises ValueError
    ],
)
def test_check_cem_strength_class(cem_class, expected):
    """Test check_cem_strength_class function."""
    try:
        assert _concrete_material_properties._check_cem_strength_class(cem_class) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            _concrete_material_properties._check_cem_strength_class(cem_class) == expected
            assert str(exc_info.value).startswith(
                    "Unknown cem_class used.")


@pytest.mark.parametrize(
    'cem_class, fcm, agg_type, expected',
    [
        ("32.5 R", 60, "basalt", {"S": 0.25, "alpha_e": 1.2}),
        ("32.5 R", 60, "quartzite", {"S": 0.25, "alpha_e": 1.0}),
        ("32.5 R", 60, "limestone", {"S": 0.25, "alpha_e": 0.9}),
        ("32.5 R", 60, "sandstone", {"S": 0.25, "alpha_e": 0.7}),
        ("32.5 R", 61, "basalt", {"S": 0.2, "alpha_e": 1.2}),
        ("32.5 R", 61, "quartzite", {"S": 0.2, "alpha_e": 1.0}),
        ("32.5 R", 61, "limestone", {"S": 0.2, "alpha_e": 0.9}),
        ("32.5 R", 61, "sandstone", {"S": 0.2, "alpha_e": 0.7}),
        ("42.5 R", 60, "basalt", {"S": 0.2, "alpha_e": 1.2}),
        ("42.5 R", 60, "Quartzite", {"S": 0.2, "alpha_e": 1.0}),
        ("42.5 R", 60, "LimeStone", {"S": 0.2, "alpha_e": 0.9}),
        ("42.5 R", 60, "sandstone", {"S": 0.2, "alpha_e": 0.7}),
        ("32.5 N", 60, "basalt", {"S": 0.38, "alpha_e": 1.2}),
        ("32.5 N", 60, "quartzite", {"S": 0.38, "alpha_e": 1.0}),
        ("32.5 N", 60, "limestone", {"S": 0.38, "alpha_e": 0.9}),
        ("32.5 N", 60, "sandstone", {"S": 0.38, "alpha_e": 0.7}),
        ("32.5 N", 61, "basalt", {"S": 0.2, "alpha_e": 1.2}),
        ("32.5 N", 61, "quartzite", {"S": 0.2, "alpha_e": 1.0}),
        ("32.5 N", 61, "limestone", {"S": 0.2, "alpha_e": 0.9}),
        ("32.5 N", 61, "sandstone", {"S": 0.2, "alpha_e": 0.7})
    ],
)
def test_get_Ecmod_coeffs(cem_class, fcm, agg_type, expected):
    """Test get_Ecmod_coeffs function."""
    assert _concrete_material_properties._get_Ecmod_coeffs(cem_class, fcm, agg_type) == expected


@pytest.mark.parametrize(
    'fc, TABULAR_VALUES, fc_value_type, expected',
    [
        (28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            "mean", 36364.06146
        ),
        (20, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            "char", 36364.06146
        ),
        (28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.0},
            "mean", 30303.38455
        ),
        (28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 0.9},
            "mean", 27273.04609
        ),
        (28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 0.7},
            "mean", 21212.36918
        )
    ],
)
def test_calc_E28(fc, TABULAR_VALUES, fc_value_type, expected):
    """Test _calc_E28 function."""
    assert np.isclose(
        _concrete_material_properties._calc_E28(fc, TABULAR_VALUES, fc_value_type),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'time, TABULAR_VALUES, expected',
    [
        (1, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            0.342022
        ),
        (7, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            0.778801
        ),
        (28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            1.
        ),
        (60, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            1.08244
        ),
    ]
)
def test_calc_beta_cc(time, TABULAR_VALUES, expected):
    """Test _calc_beta_cc function."""
    assert np.isclose(
        _concrete_material_properties._calc_beta_cc(time, TABULAR_VALUES),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'beta_cc, expected',
    [
        (0.342022, .584827),
        (0.778801, .882497),
        (1, 1),
        (1.08244, 1.04040)
    ],
)
def test_calc_beta_e(beta_cc, expected):
    """Test _calc_beta_e function."""
    assert np.isclose(
        _concrete_material_properties._calc_beta_e(beta_cc),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'time, fc, TABULAR_VALUES, fc_value_type, expected',
    [
        (1, 28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            "mean", 21266.68497
        ),
        (7, 28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            "mean", 32091.17160
        ),
        (7, 20, {"alpha": 0, "S": 0.25, "alpha_as": 700,
                 "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
         "char", 32091.17160
         ),
        (28, 28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            "mean", 36364.06146
        ),
        (60, 28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            "mean", 37833.16954
        )
    ],
)
def test_calc_E(time, fc, TABULAR_VALUES, fc_value_type, expected):
    """Test _calc_E function."""
    assert np.isclose(
        _concrete_material_properties._calc_E(time, fc, TABULAR_VALUES, fc_value_type),
        expected, rtol=1e-5
    )
