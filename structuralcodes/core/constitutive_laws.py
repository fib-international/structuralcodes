"""Collection of constitutive laws"""
import typing as t
from .base import ConstitutiveLaw


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

    def getStress(self, eps: float) -> float:
        """Return stress given strain"""
        return self._E * eps

    def getTangent(self, eps: float) -> float:
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
        self._Eh = Eh
        self._eps_su = eps_su
        self._eps_sy = fy / E

    def getStress(self, eps: float) -> float:
        """Return the stress given strain"""
        tol = 1.0e-6
        if self._eps_su is not None:
            if abs(eps) - self._eps_su > tol:
                return 0
        sign = 1 if eps >= 0 else -1
        sig = sign * min(self._fy, abs(self._E * eps))
        # TO-DO: add hardening
        return sig

    def getTangent(self, eps: float) -> float:
        """Return the tanet for given strain"""
        tol = 1.0e-6
        if abs(eps) - self._eps_sy > tol:
            return 0.0
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

    def getStress(self, eps: float) -> float:
        """Return the stress given strain"""
        sigma = 0.0
        if eps < 0.0:
            if eps >= self._eps_0:
                # Parabolic branch
                sigma = self._fc * (1 - (1 - (eps / self._eps_0)) ** self._n)
            elif eps >= self._eps_u:
                sigma = self._fc
        return sigma

    def getTangent(self, eps: float) -> float:
        """Return the tangent given strain"""
        # this function is still TO DO
        return 0.0
