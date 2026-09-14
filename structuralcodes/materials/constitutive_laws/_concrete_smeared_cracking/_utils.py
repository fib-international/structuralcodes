"""Utility functions for concrete smeared cracking."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike


def calculate_principal_strains(
    strain: ArrayLike,
) -> t.Tuple[ArrayLike, float]:
    """Calculate the principal strains in 2D."""
    eps = np.atleast_1d(strain)
    # Use arctan2 to obtain the angle in the correct quadrant
    phi = 0.5 * np.arctan2(eps[2], eps[0] - eps[1])
    eps_p = np.array(
        [
            (eps[0] + eps[1]) / 2
            + (((eps[0] - eps[1]) / 2) ** 2 + 0.25 * eps[2] ** 2) ** 0.5,
            (eps[0] + eps[1]) / 2
            - (((eps[0] - eps[1]) / 2) ** 2 + 0.25 * eps[2] ** 2) ** 0.5,
        ]
    )

    return eps_p, phi


def establish_strain_transformation_matrix(phi: float) -> ArrayLike:
    """Establish the strain transformation matrix for a given angle."""
    c = np.cos(phi)
    s = np.sin(phi)
    return np.array(
        [
            [c * c, s * s, s * c],
            [s * s, c * c, -s * c],
            [-2 * s * c, 2 * s * c, c * c - s * s],
        ]
    )
