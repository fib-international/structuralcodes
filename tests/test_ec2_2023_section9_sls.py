"""Tests for the _section9_sls module"""

import math

import pytest

from structuralcodes.codes.ec2_2023 import _section9_sls


@pytest.mark.parametrize(
    'test_input1, test_input2, expected',
    [
        (33, 1.5, 13929.86),
        (33, 2, 11608.22),
        (33, 2.5, 9949.9),
        (43, 1.5, 14829.57),
        (43, 2, 12357.98),
        (43, 2.5, 10592.55),
        (53, 1.5, 15507.43),
        (53, 2, 12922.86),
        (53, 2.5, 11076.74),
        (68, 2, 13540.73),
        (68, 2, 13540.73),
        (68, 2, 13540.73),
        (68, 2, 13540.73),
    ],
)
def test_Ec_eff(test_input1, test_input2, expected):
    """Test the Ec_eff function."""
    assert math.isclose(
        _section9_sls.Ec_eff(test_input1, test_input2), expected, rel_tol=0.01
    )


@pytest.mark.parametrize(
    'test_input1, test_input2, test_input3, test_input4, test_input5,'
    + 'expected1, expected2',
    [
        (-1000, 1, 0.15, 2.565, 500, 0, 0),
        (-500, 1, 0.15, 2.565, 500, 0, 0),
        (0, 1, 0.15, 2.565, 500, 1.231, 0),
        (500, 1, 0.15, 2.565, 500, 3.078, 3.078),
        (3000, 1, 0.15, 2.565, 500, 3.078, 3.078),
        (-1000, 1, 0.5, 2.565, 400, 0, 0),
        (-500, 1, 0.5, 2.565, 400, 0.611, 0),
        (0, 1, 0.5, 2.565, 400, 4.361, 0),
        (500, 1, 0.5, 2.565, 400, 8.111, 4.389),
        (1000, 1, 0.5, 2.565, 400, 10.901, 10.901),
        (-1000, 1, 0.5, 3.795, 500, 0, 0),
        (-500, 1, 0.8, 3.795, 500, 3.072, 0),
        (0, 1, 0.8, 3.795, 500, 6.072, 0),
        (500, 1, 0.8, 3.795, 500, 9.072, 0.928),
        (1000, 1, 0.8, 3.795, 500, 12.072, 7.928),
    ],
)
def test_As_min_y(
    test_input1,
    test_input2,
    test_input3,
    test_input4,
    test_input5,
    expected1,
    expected2,
):
    """Test the As_min_y function."""
    (As_min_y1, As_min_y2) = _section9_sls.As_min_y(
        test_input1, test_input2, test_input3, test_input4, test_input5
    )

    assert math.isclose(As_min_y1, expected1, rel_tol=0.05)
    assert math.isclose(As_min_y2, expected2, rel_tol=0.05)


@pytest.mark.parametrize(
    'test_input1, test_input2, expected',
    [
        (0.15, 0.15, 0.8),
        (0.15, 0.25, 0.8),
        (0.25, 0.15, 0.8),
        (0.25, 0.25, 0.8),
        (0.3, 0.3, 0.8),
        (0.4, 0.4, 0.74),
        (0.5, 0.5, 0.68),
        (0.6, 0.6, 0.62),
        (0.7, 0.7, 0.56),
        (0.8, 0.8, 0.5),
        (0.9, 0.9, 0.5),
        (1, 1, 0.5),
        (1.1, 1.1, 0.5),
    ],
)
def test_kh(test_input1, test_input2, expected):
    """Test the kh function."""
    assert math.isclose(
        _section9_sls.kh(test_input1, test_input2), expected, rel_tol=0.05
    )


@pytest.mark.parametrize(
    'test_input1, test_input2, test_input3, expected',
    [
        (200, 30, 45, 1.36),
        (500, 120, 45, 1.134),
        (1000, 200, 45, 1.06),
        (200, 30, 50, 1.417),
        (500, 120, 55, 1.169),
        (1000, 200, 82, 1.114),
        (300, 30, 58, 1.274),
        (450, 30, 32, 1.082),
        (650, 110, 45, 1.091),
    ],
)
def test_k_1_r(test_input1, test_input2, test_input3, expected):
    """Test the k1/r function."""
    assert math.isclose(
        _section9_sls.k_1_r(test_input1, test_input2, test_input3),
        expected,
        rel_tol=0.01,
    )


@pytest.mark.parametrize(
    'test_input1, test_input2, test_input3, test_input4, test_input5,'
    + ' test_input6, expected',
    [
        (100, 0.4, 4.072, 0.005, 5.365, 200000, 0.0003),
        (100, 0.6, 4.072, 0.005, 5.365, 200000, 0.0002),
        (150, 0.4, 2.565, 0.01, 6.354, 200000, 0.00045),
        (150, 0.6, 2.565, 0.01, 6.354, 200000, 0.0003),
        (200, 0.4, 2.896, 0.015, 6.091, 200000, 0.0006),
        (200, 0.6, 2.896, 0.015, 6.091, 200000, 0.0004),
        (250, 0.4, 3.21, 0.02, 5.869, 200000, 0.000891),
        (250, 0.6, 3.21, 0.02, 5.869, 200000, 0.000712),
        (300, 0.4, 3.795, 0.025, 5.512, 200000, 0.001155),
        (300, 0.6, 3.795, 0.025, 5.512, 200000, 0.000982),
    ],
)
def test_eps_sm_eps_cm(
    test_input1,
    test_input2,
    test_input3,
    test_input4,
    test_input5,
    test_input6,
    expected,
):
    """Test the k1/r function."""
    assert math.isclose(
        _section9_sls.epssm_epscm(
            test_input1,
            test_input2,
            test_input3,
            test_input4,
            test_input5,
            test_input6,
        ),
        expected,
        rel_tol=0.01,
    )


@pytest.mark.parametrize(
    'test_input1, test_input2, test_input3, expected',
    [
        (1000, 500, 175, 0.825),
        (1000, 250, 175, 0.883),
        (500, 250, 250, 0.5),
        (150, 75, 80, 0.5),
        (150, 0.1, 80, 0.733),
    ],
)
def test_kfl(test_input1, test_input2, test_input3, expected):
    """Test the kfl function."""
    assert math.isclose(
        _section9_sls.kfl(test_input1, test_input2, test_input3),
        expected,
        rel_tol=0.01,
    )


@pytest.mark.parametrize(
    'test_input1, test_input2, test_input3, test_input4, test_input5,'
    + ' test_input6,test_input7, test_input8, expected',
    [
        (25, 0.5, 0.9, 16, 0.015, 1.3, 1000, 15, 104.167),
        (40, 0.5, 0.9, 16, 0.015, 1.3, 1000, 15, 126.667),
        (40, 0.7, 0.9, 25, 0.015, 1.2, 1000, 15, 205.833),
        (40, 0.7, 0.9, 25, 0.025, 1.2, 1000, 15, 147.5),
        (40, 1, 0.9, 25, 0.025, 1.2, 1000, 15, 185),
        (40, 1, 0.9, 25, 0.025, 1.2, 150, 3, 159.25),
        (40, 1, 0.9, 25, 0.025, 1.3, 150, 3, 147),
        (80, 0.65, 1.2, 20, 0.025, 1.3, 250, 3, 206.667),
        (80, 0.65, 1.2, 20, 0.01, 1.3, 250, 3, 247),
    ],
)
def test_srm_cal(
    test_input1,
    test_input2,
    test_input3,
    test_input4,
    test_input5,
    test_input6,
    test_input7,
    test_input8,
    expected,
):
    """Test the srm_cal function."""
    assert math.isclose(
        _section9_sls.srm_cal(
            test_input1,
            test_input2,
            test_input3,
            test_input4,
            test_input5,
            test_input6,
            test_input7,
            test_input8,
        ),
        expected,
        rel_tol=0.01,
    )


@pytest.mark.parametrize(
    'test_input1, test_input2, test_input3, test_input4, test_input5,'
    + ' test_input6, test_input7, test_input8, test_input9, test_input10,'
    + ' test_input11, test_input12, test_input13, test_input14,'
    + ' expected1, expected2, expected3, expected4',
    [
        (
            1.3,
            1000,
            500,
            175,
            45,
            0.9,
            16,
            0.0153,
            110,
            253,
            0.4,
            2.56,
            6.35,
            200000,
            0.218,
            1.063,
            175.343,
            0.000898,
        ),
        (
            1.3,
            500,
            150,
            55,
            80,
            0.9,
            25,
            0.035,
            0.225,
            320,
            0.6,
            3.56,
            7,
            200000,
            0.394,
            1.227,
            202.27,
            0.00122,
        ),
        (
            1.2,
            500,
            250,
            55,
            80,
            0.9,
            25,
            0.035,
            0.225,
            320,
            0.6,
            3.56,
            7,
            200000,
            0.358,
            1.227,
            199.464,
            0.00122,
        ),
        (
            1.2,
            500,
            250,
            55,
            80,
            0.9,
            12,
            0.012,
            0.225,
            320,
            0.6,
            3.56,
            7,
            200000,
            0.215,
            1.208,
            231.25,
            0.00064,
        ),
        (
            1.3,
            300,
            150,
            150,
            35,
            0.9,
            16,
            0.0085382,
            0.05,
            250,
            0.4,
            2.89,
            6.09,
            200000,
            0.193,
            1.167,
            169.621,
            0.00075,
        ),
        (
            1.3,
            300,
            150,
            150,
            35,
            1.2,
            16,
            0.0085382,
            0.05,
            250,
            0.4,
            2.89,
            6.09,
            200000,
            0.237,
            1.167,
            208.661,
            0.00075,
        ),
        (
            1.3,
            300,
            150,
            150,
            45,
            1.2,
            16,
            0.0085382,
            0.075,
            250,
            0.4,
            2.89,
            6.09,
            200000,
            0.265,
            1.215,
            223.661,
            0.00075,
        ),
    ],
)
def test_wk_cal(
    test_input1,
    test_input2,
    test_input3,
    test_input4,
    test_input5,
    test_input6,
    test_input7,
    test_input8,
    test_input9,
    test_input10,
    test_input11,
    test_input12,
    test_input13,
    test_input14,
    expected1,
    expected2,
    expected3,
    expected4,
):
    """Test the wk_cal function."""
    (wk_cal_, k_1_r_, srm_cal_, epssm_epscm_) = _section9_sls.wk_cal(
        test_input1,
        test_input2,
        test_input3,
        test_input4,
        test_input5,
        test_input6,
        test_input7,
        test_input8,
        test_input9,
        test_input10,
        test_input11,
        test_input12,
        test_input13,
        test_input14,
    )

    assert math.isclose(wk_cal_, expected1, rel_tol=0.01)
    assert math.isclose(k_1_r_, expected2, rel_tol=0.01)
    assert math.isclose(srm_cal_, expected3, rel_tol=0.01)
    assert math.isclose(epssm_epscm_, expected4, rel_tol=0.01)


@pytest.mark.parametrize(
    'test_input1, test_input2,test_input3,test_input4,'
    + 'test_input5, test_input6, test_input7, test_input8,'
    + 'test_input9, expected',
    [
        (9.42, 5.32, 25, 2, 1, 0.62, 0.57, 49.1, 651.87, 24.45),
        (0.71, 0, 25, 2, 1, 0.3, 0.25, 49.1, 47.65, 0.82),
        (4, 1.85, 25, 2, 1, 0.55, 0.5, 49.1, 571.82, 9.1),
        (2.83, 1.77, 35, 2.5, 1, 0.62, 0.57, 11.2, 148.7, 4.6),
        (0.16, 0, 35, 2.5, 1, 0.2, 0.15, 11.2, 6.52, 0.16),
        (1.38, 0.75, 35, 2.5, 1, 1.05, 1, 11.2, 260.87, 2.13),
        (4.88, 3.06, 45, 2, 1, 0.62, 0.57, 31.4159265358979, 417.09, 17.54),
        (0.27, 0, 45, 2, 1, 0.2, 0.15, 31.4159265358979, 18.29, 0.27),
        (2.38, 1.29, 45, 2, 1, 1.05, 1, 31.4159265358979, 731.74, 8.15),
        (6.3, 3.39, 25, 2.5, 1, 0.62, 0.57, 20.1, 266.86, 23.54),
        (0.35, 0, 25, 2.5, 1, 0.2, 0.15, 20.1, 11.7, 0.35),
        (3.07, 1.43, 25, 2.5, 1, 1.05, 1, 20.1, 468.17, 4.51),
    ],
)
def test_delta_simpl(
    test_input1,
    test_input2,
    test_input3,
    test_input4,
    test_input5,
    test_input6,
    test_input7,
    test_input8,
    test_input9,
    expected,
):
    """Test the deslta_simpl function."""
    assert math.isclose(
        _section9_sls.delta_simpl(
            test_input1,
            test_input2,
            test_input3,
            test_input4,
            test_input5,
            test_input6,
            test_input7,
            test_input8,
            test_input9,
        ),
        expected,
        rel_tol=0.01,
    )
