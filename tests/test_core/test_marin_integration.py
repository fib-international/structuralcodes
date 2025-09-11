"""Tests for the Marin integration."""

import math

import numpy as np

from structuralcodes.core._marin_integration import marin_integration


def test_area_moments_rectangular_section():
    """Test calculating area moments of a rectangular section."""
    # Arrange
    width = 200
    height = 500

    area = width * height
    I_y = 1 / 12 * width * height**3
    I_z = 1 / 12 * height * width**3

    y = np.array([-width / 2, width / 2, width / 2, -width / 2, -width / 2])
    z = np.array(
        [-height / 2, -height / 2, height / 2, height / 2, -height / 2]
    )

    # Act
    area_marin = marin_integration(y, z, 0, 0)
    I_xy_marin = marin_integration(y, z, 1, 1)
    I_y_marin = marin_integration(y, z, 0, 2)
    I_z_marin = marin_integration(y, z, 2, 0)

    # Assert
    assert math.isclose(area_marin, area)
    assert math.isclose(I_xy_marin, 0.0)
    assert math.isclose(I_y_marin, I_y)
    assert math.isclose(I_z_marin, I_z)


def test_area_moments_hollow_rectangle():
    """Test calculating area moments of a hollow rectangular section."""
    # Arrange
    width = 1000
    height = 500
    width_opening = 500
    height_opening = 250

    area = width * height - width_opening * height_opening
    I_y = 1 / 12 * (width * height**3 - width_opening * height_opening**3)
    I_z = 1 / 12 * (width**3 * height - width_opening**3 * height_opening)

    y = np.array(
        [
            -width / 2,
            width / 2,
            width / 2,
            -width / 2,
            -width / 2,
            -width_opening / 2,
            -width_opening / 2,
            width_opening / 2,
            width_opening / 2,
            -width_opening / 2,
            -width / 2,
        ]
    )
    z = np.array(
        [
            -height / 2,
            -height / 2,
            height / 2,
            height / 2,
            -height / 2,
            -height_opening / 2,
            height_opening / 2,
            height_opening / 2,
            -height_opening / 2,
            -height_opening / 2,
            -height / 2,
        ]
    )

    # Act
    area_marin = marin_integration(y, z, 0, 0)
    I_y_marin = marin_integration(y, z, 0, 2)
    I_z_marin = marin_integration(y, z, 2, 0)

    # Assert
    assert math.isclose(area_marin, area)
    assert math.isclose(I_y_marin, I_y)
    assert math.isclose(I_z_marin, I_z)


def test_area_moments_triangle():
    """Test calculating area moments of a triangle with base b, height h, and
    top vertex displacement a. Axis colinear with the base.
    """
    # Arrange
    b = 500
    h = 600
    a = 200

    area = 0.5 * b * h
    I_y = 1 / 12 * b * h**3
    I_z = 1 / 12 * (b**3 * h + b**2 * h * a + b * h * a**2)

    y = np.array([0.0, b, a, 0.0])
    z = np.array([0.0, 0.0, h, 0.0])

    # Act
    area_marin = marin_integration(y, z, 0, 0)
    I_y_marin = marin_integration(y, z, 0, 2)
    I_z_marin = marin_integration(y, z, 2, 0)

    # Assert
    assert math.isclose(area_marin, area)
    assert math.isclose(I_y_marin, I_y)
    assert math.isclose(I_z_marin, I_z)
