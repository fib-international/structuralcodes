"""Collection of some standard constitutive laws"""
import typing as t
import numpy as np
from ..core.base import ConstitutiveLaw


class Elastic(ConstitutiveLaw):
    """Class for elastic Constitutive Law"""

    __materials__: t.Tuple[str] = (
        'concrete',
        'steel',
        'rebars',
    )

    def __init__(self, E: float, name: t.Optional[str] = None) -> None:
        """Initialize an Elastic Material"""
        name = name if name is not None else "ElasticLaw"
        super().__init__(name=name)
        self._E = E

    def get_stress(self, eps: float) -> float:
        """Return stress given strain"""
        eps = np.asarray(eps)
        return self._E * eps

    def get_tangent(self, eps: float) -> float:
        """Return the tangent"""
        return self._E


class ElasticPlastic(ConstitutiveLaw):
    """Class for elastic-plastic Constitutive Law"""

    __materials__: t.Tuple[str] = (
        'steel',
        'rebars',
    )

    def __init__(
        self,
        E: float,
        fy: float,
        Eh: t.Optional[float] = None,
        eps_su: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize an Elastic-Plastic Material"""
        name = name if name is not None else "ElasticPlasticLaw"
        super().__init__(name=name)
        self._E = E
        self._fy = fy
        if Eh is not None:
            self._Eh = Eh
        else:
            self._Eh = 0.0
        self._eps_su = eps_su
        self._eps_sy = fy / E

    def get_stress(self, eps: float) -> float:
        """Return the stress given strain"""
        eps = np.atleast_1d(np.asarray(eps))
        sig = self._E * eps
        delta_sig = self._fy * (1 - self._Eh / self._E)
        sig[sig < -self._fy] = eps[sig < -self._fy] * self._Eh - delta_sig
        sig[sig > self._fy] = eps[sig > self._fy] * self._Eh + delta_sig
        if self._eps_su is not None:
            sig[eps > self._eps_su] = 0
            sig[eps < -self._eps_su] = 0  # pylint: disable=E1130
        return sig

    def get_tangent(self, eps: float) -> float:
        """Return the tanet for given strain"""
        tol = 1.0e-6
        if abs(eps) - self._eps_sy > tol:
            return self._Eh
        return self._E


class ParabolaRectangle(ConstitutiveLaw):
    """Class for parabola rectangle constitutive law"""

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_0: float = -0.002,
        eps_u: float = -0.0035,
        n: float = 2.0,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize an Elastic-Plastic Material"""
        name = name if name is not None else "ParabolaRectangleLaw"
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._eps_0 = -abs(eps_0)
        self._eps_u = -abs(eps_u)
        self._n = n

    def get_stress(self, eps: float) -> float:
        """Return the stress given strain"""
        eps = np.atleast_1d(np.asarray(eps))
        # Parabolic branch
        sig = self._fc * (1 - (1 - (eps / self._eps_0)) ** self._n)
        # Rectangle branch
        sig[eps < self._eps_0] = self._fc
        # Zero elsewhere
        sig[eps < self._eps_u] = 0
        sig[eps > 0] = 0
        return sig

    def get_tangent(self, eps: float) -> float:
        """Return the tangent given strain"""
        # this function is still TO DO
        return 0.0


class UserDefined(ConstitutiveLaw):
    """Class for a user defined constitutive law
    The curve is defined with positive and optionally negative
    values. After the last value, the stress is mantained constant"""

    __materials__: t.Tuple[str] = ('concrete', 'steel', 'rebars')

    def __init__(
        self, x: t.List[float], y: t.List[float], name: t.Optional[str] = None
    ) -> None:
        """Initialize a UserDefined constitutive law"""
        name = name if name is not None else "UserDefinedLaw"
        super().__init__(name=name)
        x = np.atleast_1d(np.asarray(x))
        y = np.atleast_1d(np.asarray(y))
        if len(x) != len(y):
            raise ValueError("The two arrays should have the same length")
        if not np.any(x < 0):
            # User provided only positive part, reflect in negative
            self._x = np.concatenate((-np.flip(x), x))
            self._y = np.concatenate((-np.flip(y), y))
        else:
            # User gave both positive and negative parts
            self._x = x
            self._y = y

    def get_stress(self, eps: float) -> float:
        """Return the stress given strain"""
        eps = np.atleast_1d(np.asarray(eps))
        sig = np.interp(eps, self._x, self._y, left=0, right=0)
        return sig

    def get_tangent(self, eps: float) -> float:
        """Return the tangent given strain"""
        # this function is still TO DO
        return 0.0
