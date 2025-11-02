"""User defined constitutive law."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class UserDefined(ConstitutiveLaw):
    """Class for a user defined constitutive law.

    The curve is defined with positive and optionally negative values. After
    the last value, the stress can go to zero to simulate failure (default), or
    be maintained constante, or the last tanget or secant values may be
    maintained indefinetely. The flag parameter controls this behavior.
    """

    __materials__: t.Tuple[str] = ('concrete', 'steel', 'rebars')

    def __init__(
        self,
        x: ArrayLike,
        y: ArrayLike,
        name: t.Optional[str] = None,
        eps_u: t.Optional[t.Union[float, t.Tuple[float, float]]] = None,
        flag: int = 0,
    ) -> None:
        """Initialize a UserDefined constitutive law.

        Arguments:
            x (ArrayLike): Data for strains.
            y (ArrayLike): Data for stresses. Must be of same length as x.

        Keyword Arguments:
            name (Optional, str): A name for the constitutive law.
            eps_u (float or (float, float), optional): Ultimate strain.
                If a single value is provided the same is adopted for both
                negative and positive strains. If a tuple is provided, it
                should be given as (negative, positive). Default value = None.
            flag (Optional): A flag specifying the behavior after the last
                point. Admissible values: 0 (default): stress drops to zero
                after ultimate strain, 1: stress is mantained constant, 2:
                last tangent is used, 3: last secant is used.
        """
        name = name if name is not None else 'UserDefinedLaw'
        super().__init__(name=name)
        x = np.atleast_1d(np.asarray(x))
        y = np.atleast_1d(np.asarray(y))
        if len(x) != len(y):
            raise ValueError('The two arrays should have the same length')
        if not np.any(x < 0):
            # User provided only positive part, reflect in negative
            self._x = np.concatenate((-np.flip(x)[:-1], x))
            self._y = np.concatenate((-np.flip(y)[:-1], y))
        else:
            # User gave both positive and negative parts
            self._x = x
            self._y = y
        # Define what happens after last strain
        if flag not in (0, 1, 2, 3):
            raise ValueError('Flag can assume values 0, 1, 2 or 3.')
        self._ultimate_strain_p = self._x[-1]
        self._ultimate_strain_n = self._x[0]
        if flag in (1, 2, 3):
            x = np.insert(self._x, 0, self._x[0] * 100)
            x = np.append(x, self._x[-1] * 100)
            if flag == 1:
                y = np.insert(self._y, 0, self._y[0])
                y = np.append(y, self._y[-1])
            elif flag == 2:
                tangent_p = (self._y[-1] - self._y[-2]) / (
                    self._x[-1] - self._x[-2]
                )
                tangent_n = (self._y[1] - self._y[0]) / (
                    self._x[1] - self._x[0]
                )
                y = np.insert(
                    self._y, 0, (x[0] - x[1]) * tangent_n + self._y[0]
                )
                y = np.append(y, (x[-1] - x[-2]) * tangent_p + self._y[-1])
            elif flag == 3:
                secant_p = self._y[-1] / self._x[-1]
                secant_n = self._y[0] / self._x[0]
                y = np.insert(
                    self._y, 0, (x[0] - x[1]) * secant_n + self._y[0]
                )
                y = np.append(y, (x[-1] - x[-2]) * secant_p + self._y[-1])
            self._x = x
            self._y = y

        # Compute slope of each segment
        self._slopes = np.diff(self._y) / np.diff(self._x)

        # Set ultimate strains if needed
        if eps_u is not None:
            self._set_ultimate_strain(eps_u=eps_u)

    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the stress given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        return np.interp(eps, self._x, self._y, left=0, right=0)

    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the tangent given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)

        # Find the segment index for each x value
        indices = np.searchsorted(self._x, eps) - 1

        # Check that indices are within vlaid range
        indices = np.clip(indices, 0, len(self._slopes) - 1)

        # Get the corresponding slopes
        tangent = self._slopes[indices]

        # Elsewhere tangent is zero
        if np.isscalar(eps):
            if eps < self._x[0] or eps > self._x[-1]:
                tangent = 0
        else:
            tangent[eps < self._x[0]] = 0.0
            tangent[eps > self._x[-1]] = 0.0

        return tangent

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
        strains = []
        coeff = []
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # understand in which branch are we
            strain[0] = self.preprocess_strains_with_limits(strain[0])
            found = False
            for i in range(len(self._x) - 1):
                if self._x[i] <= strain[0] and self._x[i + 1] >= strain[0]:
                    strains = None
                    stiffness = (self._y[i + 1] - self._y[i]) / (
                        self._x[i + 1] - self._x[i]
                    )
                    a0 = stiffness * (strain[0] - self._x[i]) + self._y[i]
                    a1 = stiffness * strain[1]
                    coeff.append((a0, a1))
                    found = True
                    break
            if not found:
                strains = None
                coeff.append((0.0,))
        else:
            for i in range(len(self._x) - 1):
                # For each branch of the linear piecewise function
                stiffness = (self._y[i + 1] - self._y[i]) / (
                    self._x[i + 1] - self._x[i]
                )
                strains.append((self._x[i], self._x[i + 1]))
                a0 = stiffness * (strain[0] - self._x[i]) + self._y[i]
                a1 = stiffness * strain[1]
                coeff.append((a0, a1))

        return strains, coeff

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
        strains = []
        coeff = []
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # understand in which branch are we
            strain[0] = self.preprocess_strains_with_limits(strain[0])
            found = False
            for i in range(len(self._x) - 1):
                if self._x[i] <= strain[0] and self._x[i + 1] >= strain[0]:
                    strains = None
                    stiffness = (self._y[i + 1] - self._y[i]) / (
                        self._x[i + 1] - self._x[i]
                    )
                    coeff.append((stiffness,))
                    found = True
                    break
            if not found:
                strains = None
                coeff.append((0.0,))
        else:
            for i in range(len(self._x) - 1):
                # For each branch of the linear piecewise function
                stiffness = (self._y[i + 1] - self._y[i]) / (
                    self._x[i + 1] - self._x[i]
                )
                strains.append((self._x[i], self._x[i + 1]))
                coeff.append((stiffness,))

        return strains, coeff

    def get_ultimate_strain(self, **kwargs) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        del kwargs
        return (self._ultimate_strain_n, self._ultimate_strain_p)

    def _set_ultimate_strain(
        self, eps_u=t.Union[float, t.Tuple[float, float]]
    ) -> None:
        """Set ultimate strains for Elastic Material if needed.

        Arguments:
            eps_u (float or (float, float)): Defining ultimate strain. If a
                single value is provided the same is adopted for both negative
                and positive strains. If a tuple is provided, it should be
                given as (negative, positive).
        """
        if isinstance(eps_u, float):
            self._ultimate_strain_p = abs(eps_u)
            self._ultimate_strain_n = -abs(eps_u)
        elif isinstance(eps_u, tuple):
            if len(eps_u) < 2:
                raise ValueError(
                    'Two values need to be provided when setting the tuple'
                )
            eps_u_n = eps_u[0]
            eps_u_p = eps_u[1]
            if eps_u_p < eps_u_n:
                eps_u_p, eps_u_n = eps_u_n, eps_u_p
            if eps_u_p < 0:
                raise ValueError(
                    'Positive ultimate strain should be non-negative'
                )
            if eps_u_n > 0:
                raise ValueError(
                    'Negative utimate strain should be non-positive'
                )
            self._ultimate_strain_p = eps_u_p
            self._ultimate_strain_n = eps_u_n
        else:
            raise ValueError(
                'set_ultimate_strain requires a single value or a tuple \
                with  two values'
            )
