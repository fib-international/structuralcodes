"""Parallel wrapper constitutive law."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class Parallel(ConstitutiveLaw):
    """Class for a Parallel Constitutive Law."""

    _strain_compatibility: bool = True

    __materials__: t.Tuple[str] = (
        'steel',
        'rebars',
        'concrete',
    )

    _wrapped_law: ConstitutiveLaw = None

    def __init__(
        self,
        constitutive_laws: t.List[ConstitutiveLaw],
        weights: t.Optional[t.List[float]] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a Parallel Constitutive Law.

        This constitutive law is a wrapper for other constitutive laws
        that are assumed to be working in parallel.
        The strains are equal while stresses and stiffnesses are additive.

        Arguments:
            constitutive_laws (List[ConstitutiveLaw]): List of wrapped
                constitutive laws that are working in parallel.
            weights (List[float], optional): List of weights for each
                constitutive law. If None (Default), all weights are set to 1.
        """
        name = name if name is not None else 'ParallelLaw'
        super().__init__(name=name)
        for idx, law in enumerate(constitutive_laws):
            if not isinstance(law, ConstitutiveLaw):
                raise TypeError(
                    f'Element at index {idx} is of type {type(law)}, expected',
                    ' ConstitutiveLaw',
                )
        self._wrapped_laws = constitutive_laws

        if weights is None:
            self._weights = np.ones(len(constitutive_laws))
        else:
            if len(weights) != len(constitutive_laws):
                raise ValueError(
                    'Length of weights must be equal to length of '
                    'constitutive_laws'
                )
            self._weights = np.array(weights)

    @property
    def weights(self) -> np.ndarray:
        """Return the weights of the wrapped constitutive laws."""
        return self._weights

    @property
    def wrapped_laws(self) -> t.List[ConstitutiveLaw]:
        """Return the wrapped constitutive laws."""
        return self._wrapped_laws

    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the stress given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        stress = 0 * eps  # This can be scalar or array
        for weight, law in zip(self._weights, self._wrapped_laws):
            stress += weight * law.get_stress(eps)

        return stress

    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the tangent for given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        tangent = 0 * eps  # This can be scalar or array
        for weight, law in zip(self._weights, self._wrapped_laws):
            tangent += weight * law.get_tangent(eps)
        return tangent

    def _split_marin(self, laws_marin):
        """Split the marin data into ranges and coeffs."""
        # Compute combined intervals

        # Collect all possible breakpoints in a set to avoid duplicates
        def _collect_breakpoints(laws_marin):
            breakpoints = set()

            for ranges, _ in laws_marin:
                if ranges is not None:
                    for r in ranges:
                        eps0, eps1 = r
                        breakpoints.add(eps0)
                        breakpoints.add(eps1)

            return sorted(breakpoints)  # Sort them

        breakpoints = _collect_breakpoints(laws_marin)

        # Compute refined ranges
        refined_ranges = [
            (breakpoints[i], breakpoints[i + 1])
            for i in range(len(breakpoints) - 1)
        ]

        # Compute a map of coefficients for each law and each refined range
        def _map_coeffs(laws_marin, refined_ranges):
            law_maps = []
            for ranges, coeffs in laws_marin:
                if ranges is None:
                    # The law is valid for the whole range
                    if len(refined_ranges) == 0:
                        # No ranges at all, means no contribution
                        law_map = [tuple(v for v in coeffs[0])]
                    else:
                        law_map = [tuple(v for v in coeffs[0])] * len(
                            refined_ranges
                        )
                else:
                    # The law is valid for specific ranges

                    # The ranges returned by __marin__ could be not sorted
                    zipped = sorted(
                        zip(ranges, coeffs), key=lambda rc: rc[0][0]
                    )
                    sorted_ranges, sorted_coeffs = zip(*zipped)

                    law_map = []
                    j = 0
                    for eps0, eps1 in refined_ranges:
                        # Advance j till finding the matching range
                        while (
                            j < len(sorted_ranges)
                            and eps1 > sorted_ranges[j][1]
                        ):
                            j += 1
                        if (
                            j < len(sorted_ranges)
                            and eps0 >= sorted_ranges[j][0]
                            and eps1 <= sorted_ranges[j][1]
                        ):
                            law_map.append(sorted_coeffs[j])
                        else:
                            law_map.append(())  # No contribution in this range
                law_maps.append(law_map)
            return law_maps

        law_maps = _map_coeffs(laws_marin, refined_ranges)

        # Finally computed the weighted sum of coefficients for each range
        coeffs = []
        for i in range(max(len(refined_ranges), 1)):
            coeffs_law = []
            for w, law_map in zip(self.weights, law_maps):
                coeff = [w * c for c in law_map[i]] if law_map[i] else []
                coeffs_law.append(coeff)

            # Compute the max order of the polynomial for the range
            max_len = max((len(c) for c in coeffs_law), default=0)
            total = [0.0] * max_len
            for c in coeffs_law:
                for k in range(len(c)):
                    total[k] += c[k]
            coeffs.append(tuple(total))
        if len(refined_ranges) == 0:
            refined_ranges = None
        return (refined_ranges, coeffs)

    def __marin__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
        """Returns coefficients and strain limits for Marin integration in a
        simply formatted way.

        Arguments:
            strain (float, float): Tuple defining the strain profile: eps =
                strain[0] + strain[1]*y.

        Example:
            [(0, -0.002), (-0.002, -0.003)]
            [(a0, a1, a2), (a0)]
        """
        # Apply interval merging with piecewise polynomial combination
        laws_marin = []
        for law in self.wrapped_laws:
            laws_marin.append(law.__marin__(strain))

        return self._split_marin(laws_marin)

    def __marin_tangent__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
        """Returns coefficients and strain limits for Marin integration of
        tangent in a simply formatted way.

        Arguments:
            strain (float, float): Tuple defining the strain profile: eps =
                strain[0] + strain[1]*y.

        Example:
            [(0, -0.002), (-0.002, -0.003)]
            [(a0, a1, a2), (a0)]
        """
        # Apply interval merging with piecewise polynomial combination
        laws_marin = []
        for law in self.wrapped_laws:
            laws_marin.append(law.__marin_tangent__(strain))

        return self._split_marin(laws_marin)

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        ult_strain = [np.inf, -np.inf]
        for law in self._wrapped_laws:
            eps_u = law.get_ultimate_strain(yielding=yielding)
            ult_strain[0] = min(ult_strain[0], eps_u[0])
            ult_strain[1] = max(ult_strain[1], eps_u[1])
        return tuple(ult_strain)
