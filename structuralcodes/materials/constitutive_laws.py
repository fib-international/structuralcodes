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
        self._eps_su = None

    def get_stress(self, eps: ArrayLike) -> float:
        """Return stress given strain."""
        eps = np.asarray(eps)
        return self._E * eps

    def get_tangent(self, *args) -> float:
        """Return the tangent."""
        del args
        return self._E

    def __marin__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
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
        strains = None
        y_na = -strain[0] / strain[1]
        a0 = -self._E * y_na * strain[1]
        a1 = self._E * strain[1]
        coeff = [(a0, a1)]
        return strains, coeff

    def get_ultimate_strain(self, **kwargs) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        # There is no real strain limit, so set it to very large values
        # unlesse specified by the user differently
        del kwargs
        return self._eps_su or (100, -100)

    def set_ultimate_strain(
        self, eps_su=t.Union[float, t.Tuple[float, float]]
    ) -> None:
        """Set ultimate strains for Elastic Material if needed.

        Args:
            eps_su: float or (float, float) defining ultimate strain
                if a single value is provided the same is adopted for
                both positive and negative strains
        """
        if isinstance(eps_su, float):
            self._eps_su = (abs(eps_su), -abs(eps_su))
        elif isinstance(eps_su, tuple):
            if len(eps_su) < 2:
                raise ValueError(
                    'Two values need to be provided when setting the tuple'
                )
            eps_su_p = eps_su[0]
            eps_su_n = eps_su[1]
            if eps_su_p < eps_su_n:
                eps_su_p, eps_su_n = eps_su_n, eps_su_p
            if eps_su_p < 0:
                raise ValueError(
                    'Positive ultimate strain should be non-negative'
                )
            if eps_su_n > 0:
                raise ValueError(
                    'Negative utimate strain should be non-positive'
                )
            self._eps_su = (eps_su_p, eps_su_n)
        else:
            raise ValueError(
                'set_ultimate_strain requires a single value or a tuple \
                with  two values'
            )


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

    def __marin__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
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
        y_na = -strain[0] / strain[1]
        eps_sy_p, eps_sy_n = self.get_ultimate_strain(yielding=True)
        # y_yp = (strain[0] - eps_sy_p) / strain[1]
        # y_yn = (strain[0] - eps_sy_n) / strain[1]
        # Hardening part negative
        eps_su_p, eps_su_n = self.get_ultimate_strain()
        strains.append((eps_su_n, eps_sy_n))
        a0 = -self._Eh * y_na * strain[1] - self._fy * (1 - self._Eh / self._E)
        a1 = self._Eh * strain[1]
        coeff.append((a0, a1))
        # Elastic part
        strains.append((eps_sy_n, eps_sy_p))
        a0 = -self._E * y_na * strain[1]
        a1 = self._E * strain[1]
        coeff.append((a0, a1))
        # Hardening part positive
        strains.append((eps_sy_p, eps_su_p))
        a0 = -self._Eh * y_na * strain[1] + self._fy * (1 - self._Eh / self._E)
        a1 = self._Eh * strain[1]
        coeff.append((a0, a1))
        return strains, coeff

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        if yielding:
            return (self._eps_sy, -self._eps_sy)
        # If not specified eps
        if self._eps_su is None:
            return (self._eps_sy * 2, -self._eps_sy * 2)
        return (self._eps_su, -self._eps_su)


class ParabolaRectangle(ConstitutiveLaw):
    """Class for parabola rectangle constitutive law.
    The stresses and strains are assumed negative in compression and positive
    in tension.
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_0: float = -0.002,
        eps_u: float = -0.0035,
        n: float = 2.0,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a Parabola-Rectangle Material.

        Arguments:
        fc: (float) the strength of concrete in compression
        eps_0: (float) peak strain of concrete in compression
               optional, default value = -0.002
        eps_u: (float) ultimate strain of concrete in compression
               optional, default value = -0.0035
        n: (float) exponent for the pre-peak branch
           optional, default value = 2
        name: (str) a name for the constitutive law, optional
        """
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

    def __marin__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[float], t.List[float]]:
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
        y_na = -strain[0] / strain[1]
        y_0 = (self._eps_0 - strain[0]) / strain[1]
        # (strain[0] - self._eps_u) / strain[1]
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

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        if yielding:
            return (100, self._eps_0)
        return (100, self._eps_u)


class Sargin(ConstitutiveLaw):
    """Class for Sargin constitutive law.
    The stresses and strains are assumed negative in compression
    and positive in tension.
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_c1: float = -0.0023,
        eps_cu1: float = -0.0035,
        k: float = 2.04,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a Sargin Material.

        Arguments:
            fc: (float) the strength of concrete in compression
            eps_c1: (float) peak strain of concrete in compression
                optional, default value = -0.0023
            eps_u: (float) ultimate strain of concrete in compression
                optional, default value = -0.0035
            k: (float) plasticity numberm, default value = 2.04
            name: (str) a name for the constitutive law,

        Raises:
            ValueError: if k is less or equal to 0
        Notes:
            if positive values are input for fc, eps_c1 and eps_cu1
            are input, they will be assumed negative.
        """
        name = name if name is not None else 'ParabolaRectangleLaw'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._eps_c1 = -abs(eps_c1)
        self._eps_cu1 = -abs(eps_cu1)
        self._k = k

    def get_stress(self, eps: float) -> float:
        """Return the stress given the strain."""
        eps = np.atleast_1d(np.asarray(eps))
        # Polynomial branch
        eta = eps / self._eps_c1

        sig = self._fc * (self._k * eta - eta**2) / (1 + (self._k - 2) * eta)

        # Elsewhere stress is 0.0
        sig[eps < self._eps_cu1] = 0.0
        sig[eps > 0] = 0.0

        return sig

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        # polynomial branch
        eta = eps / self._eps_c1

        tangent = (
            self._fc
            / self._eps_c1
            * ((2 - self._k) * eta**2 - 2 * eta + self._k)
            / (1 + (self._k - 2) * eta) ** 2
        )
        # Elsewhere tangent is zero
        tangent[eps < self._eps_cu1] = 0.0
        tangent[eps > 0] = 0.0

        return tangent

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        if yielding:
            return (100, self._eps_c1)
        return (100, self._eps_cu1)


class UserDefined(ConstitutiveLaw):
    """Class for a user defined constitutive law
    The curve is defined with positive and optionally negative
    values. After the last value, the stress can go to zero to simulate
    failure (default), or be mantained constante, or the last tanget or
    secant values may be mantained indefinetely. The flag parameter
    controls this behavior.
    """

    __materials__: t.Tuple[str] = ('concrete', 'steel', 'rebars')

    def __init__(
        self,
        x: ArrayLike,
        y: ArrayLike,
        name: t.Optional[str] = None,
        flag: int = 0,
    ) -> None:
        """Initialize a UserDefined constitutive law.

        Arguments:
            x, y: two arrayLike objects containing data for strain and stress
                    (must be of same length)
            name: (Optional) a name for the constitutive law
            flag: (Optional) a flag specifying the behavior after the last
                point.
                Admissible values:
                    0 (default): stress drops to zero after ultimate strain
                    1: stress is mantained constant
                    2: last tangent is used
                    3: last secant is used
        """
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

    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Return the stress given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        return np.interp(eps, self._x, self._y, left=0, right=0)

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent given strain."""
        # this function is still TO DO
        raise NotImplementedError

    def __marin__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
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

    def get_ultimate_strain(self, **kwargs) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        del kwargs
        return (self._ultimate_strain_p, self._ultimate_strain_n)

    def set_ultimate_strain(
        self, eps_su=t.Union[float, t.Tuple[float, float]]
    ) -> None:
        """Set ultimate strains for Elastic Material if needed.

        Args:
            eps_su: float or (float, float) defining ultimate strain
                if a single value is provided the same is adopted for
                both positive and negative strains
        """
        if isinstance(eps_su, float):
            self._ultimate_strain_p = abs(eps_su)
            self._ultimate_strain_n = -abs(eps_su)
        elif isinstance(eps_su, tuple):
            if len(eps_su) < 2:
                raise ValueError(
                    'Two values need to be provided when setting the tuple'
                )
            eps_su_p = eps_su[0]
            eps_su_n = eps_su[1]
            if eps_su_p < eps_su_n:
                eps_su_p, eps_su_n = eps_su_n, eps_su_p
            if eps_su_p < 0:
                raise ValueError(
                    'Positive ultimate strain should be non-negative'
                )
            if eps_su_n > 0:
                raise ValueError(
                    'Negative utimate strain should be non-positive'
                )
            self._ultimate_strain_p = eps_su_p
            self._ultimate_strain_n = eps_su_n
        else:
            raise ValueError(
                'set_ultimate_strain requires a single value or a tuple \
                with  two values'
            )
