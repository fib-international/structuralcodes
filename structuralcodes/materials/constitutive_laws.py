"""Collection of some standard constitutive laws."""
import typing as t
import numpy as np
from numpy.typing import ArrayLike
from ..core.base import ConstitutiveLaw


class Elastic(ConstitutiveLaw):
    """Class for elastic Constitutive Law."""

    __materials__: t.Tuple[str] = (
        'concrete',
        'steel',
        'rebars',
    )

    def __init__(self, E: float, name: t.Optional[str] = None) -> None:
        """Initialize an Elastic Material."""
        name = name if name is not None else 'ElasticLaw'
        super().__init__(name=name)
        self._E = E

    def get_stress(self, eps: ArrayLike) -> float:
        """Return stress given strain."""
        eps = np.asarray(eps)
        return self._E * eps

    def get_tangent(self, eps: ArrayLike) -> float:
        """Return the tangent."""
        return self._E


class ElasticPlastic(ConstitutiveLaw):
    """Class for elastic-plastic Constitutive Law."""

    __materials__: t.Tuple[str] = (
        'steel',
        'rebars',
    )

    def __init__(
        self,
        E: float,
        fy: float,
        Eh: float = 0.0,
        eps_su: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize an Elastic-Plastic Material."""
        name = name if name is not None else 'ElasticPlasticLaw'
        super().__init__(name=name)
        if E > 0:
            self._E = E
        else:
            raise ValueError('Elastic modulus E must be greater than zero')
        self._fy = fy
        self._Eh = Eh
        self._eps_su = eps_su
        self._eps_sy = fy / E

    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Return the stress given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        sig = self._E * eps
        delta_sig = self._fy * (1 - self._Eh / self._E)
        sig[sig < -self._fy] = eps[sig < -self._fy] * self._Eh - delta_sig
        sig[sig > self._fy] = eps[sig > self._fy] * self._Eh + delta_sig
        if self._eps_su is not None:
            sig[eps > (self._eps_su * 1.01)] = 0
            sig[eps < (-self._eps_su * 1.01)] = 0  # pylint: disable=E1130
        return sig

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent for given strain."""
        tol = 1.0e-6
        if abs(eps) - self._eps_sy > tol:
            return self._Eh
        return self._E

    def get_ultimate_strain(self) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        return (self._eps_su, -self._eps_su)


class ParabolaRectangle(ConstitutiveLaw):
    """Class for parabola rectangle constitutive law."""

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_0: float = -0.002,
        eps_u: float = -0.0035,
        n: float = 2.0,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize an Elastic-Plastic Material."""
        name = name if name is not None else 'ParabolaRectangleLaw'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._eps_0 = -abs(eps_0)
        self._eps_u = -abs(eps_u)
        self._n = n

    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Return the stress given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        # Parabolic branch
        sig = self._fc * (1 - (1 - (eps / self._eps_0)) ** self._n)
        # Rectangle branch
        sig[eps < self._eps_0] = self._fc
        # Zero elsewhere
        sig[eps < self._eps_u] = 0
        sig[eps > 0] = 0
        return sig

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        # parabolic branch
        tangent = (
            self._n
            * self._fc
            / self._eps_0
            * (1 - (eps / self._eps_0)) ** (self._n - 1)
        )
        # Elsewhere tangent is zero
        tangent[eps < self._eps_0] = 0.0
        tangent[eps > 0] = 0.0
        return tangent

    def __marin__(self, strain):
        """Returns coefficients and strain limits for Marin
        integration in a simply formatted way.

        Args:
            strain: (float, float) tuple defining the strain
                profile: eps = strain[0] + strain[1]*y
        Returns:

        Example:
            [(0, -0.002), (-0.002, -0.003)]
            [(a0, a1, a2), (a0)]
        """
        strains = []
        coeff = []
        y_na = strain[0] / strain[1]
        y_0 = (strain[0] - self._eps_0) / strain[1]
        (strain[0] - self._eps_u) / strain[1]
        # Parabolic part
        strains.append((self._eps_0, 0))
        y0na = y_0 - y_na
        yna_y0na = y_na / y0na
        a0 = -2 * self._fc * yna_y0na * (1 + yna_y0na / 2.0)
        a1 = 2 * self._fc / y0na * (1 + yna_y0na)
        a2 = -self._fc / y0na**2
        coeff.append((a0, a1, a2))
        # Constant part
        strains.append((self._eps_u, self._eps_0))
        coeff.append((self._fc,))
        return strains, coeff

    def get_ultimate_strain(self) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        return (100, self._eps_u)


class UserDefined(ConstitutiveLaw):
    """Class for a user defined constitutive law
    The curve is defined with positive and optionally negative
    values. After the last value, the stress goes to zero simulating failure.
    """

    __materials__: t.Tuple[str] = ('concrete', 'steel', 'rebars')

    def __init__(
        self, x: ArrayLike, y: ArrayLike, name: t.Optional[str] = None
    ) -> None:
        """Initialize a UserDefined constitutive law."""
        name = name if name is not None else 'UserDefinedLaw'
        super().__init__(name=name)
        x = np.atleast_1d(np.asarray(x))
        y = np.atleast_1d(np.asarray(y))
        if len(x) != len(y):
            raise ValueError('The two arrays should have the same length')
        if not np.any(x < 0):
            # User provided only positive part, reflect in negative
            self._x = np.concatenate((-np.flip(x), x))
            self._y = np.concatenate((-np.flip(y), y))
        else:
            # User gave both positive and negative parts
            self._x = x
            self._y = y

    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Return the stress given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        return np.interp(eps, self._x, self._y, left=0, right=0)

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent given strain."""
        # this function is still TO DO
        raise NotImplementedError
