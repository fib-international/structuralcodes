"""Collection of some standard constitutive laws."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ..core.base import ConstitutiveLaw, Material


class Elastic(ConstitutiveLaw):
    """Class for elastic constitutive law."""

    __materials__: t.Tuple[str] = (
        'concrete',
        'steel',
        'rebars',
    )

    def __init__(self, E: float, name: t.Optional[str] = None) -> None:
        """Initialize an Elastic Material.

        Arguments:
            E (float): The elastic modulus.

        Keyword Arguments:
            name (str): A descriptive name for the constitutive law.
        """
        name = name if name is not None else 'ElasticLaw'
        super().__init__(name=name)
        self._E = E
        self._eps_su = None

    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Return stress given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        return self._E * eps

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent."""
        eps = np.atleast_1d(np.asarray(eps))
        return np.ones_like(eps) * self._E

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
        strains = None
        a0 = self._E * strain[0]
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

        Arguments:
            eps_su (float or (float, float)): Defining ultimate strain if a
                single value is provided the same is adopted for both positive
                and negative strains.
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
        """Initialize an Elastic-Plastic Material.

        Arguments:
            E (float): The elastic modulus.
            fy (float): The yield strength.

        Keyword Arguments:
            Eh (float): The hardening modulus.
            eps_su (float): The ultimate strain.
            name (str): A descriptive name for the constitutive law.
        """
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
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        sig = self._E * eps
        delta_sig = self._fy * (1 - self._Eh / self._E)
        sig[sig < -self._fy] = eps[sig < -self._fy] * self._Eh - delta_sig
        sig[sig > self._fy] = eps[sig > self._fy] * self._Eh + delta_sig
        if self._eps_su is not None:
            sig[eps > self._eps_su] = 0
            sig[eps < -self._eps_su] = 0  # pylint: disable=E1130
        return sig

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent for given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        tangent = np.ones_like(eps) * self._E
        tangent[eps > self._eps_sy] = self._Eh
        tangent[eps < -self._eps_sy] = self._Eh

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
        eps_sy_p, eps_sy_n = self.get_ultimate_strain(yielding=True)
        eps_su_p, eps_su_n = self.get_ultimate_strain()
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # Understand in which branch are we
            strain[0] = self.preprocess_strains_with_limits(strain[0])[0]
            if strain[0] > eps_sy_p and strain[0] <= eps_su_p:
                # We are in the Hardening part positive
                strains = None
                a0 = self._Eh * strain[0] + self._fy * (1 - self._Eh / self._E)
                a1 = self._Eh * strain[1]
                coeff.append((a0, a1))
            elif strain[0] < eps_sy_n and strain[0] >= eps_su_n:
                # We are in the Hardening part negative
                strains = None
                a0 = self._Eh * strain[0] - self._fy * (1 - self._Eh / self._E)
                a1 = self._Eh * strain[1]
                coeff.append((a0, a1))
            elif abs(strain[0]) <= self._eps_sy:
                # We are in the elastic part
                strains = None
                a0 = self._E * strain[0]
                a1 = self._E * strain[1]
                coeff.append((a0, a1))
            else:
                strains = None
                coeff.append((0.0,))
        else:
            # Hardening part negative
            strains.append((eps_su_n, eps_sy_n))
            a0 = self._Eh * strain[0] - self._fy * (1 - self._Eh / self._E)
            a1 = self._Eh * strain[1]
            coeff.append((a0, a1))
            # Elastic part
            strains.append((eps_sy_n, eps_sy_p))
            a0 = self._E * strain[0]
            a1 = self._E * strain[1]
            coeff.append((a0, a1))
            # Hardening part positive
            strains.append((eps_sy_p, eps_su_p))
            a0 = self._Eh * strain[0] + self._fy * (1 - self._Eh / self._E)
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
            fc (float): The strength of concrete in compression.

        Keyword Arguments:
            eps_0 (float): Peak strain of concrete in compression. Default
                value = -0.002.
            eps_u (float): Ultimate strain of concrete in compression. Default
                value = -0.0035.
            n (float): Exponent for the pre-peak branch. Default value = 2.
            name (str): A name for the constitutive law.
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
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        sig = np.zeros_like(eps)
        # Parabolic branch
        sig[(eps <= 0) & (eps >= self._eps_0)] = self._fc * (
            1
            - (1 - (eps[(eps <= 0) & (eps >= self._eps_0)] / self._eps_0))
            ** self._n
        )
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
        tangent = np.zeros_like(eps)
        tangent[(eps <= 0) & (eps >= self._eps_0)] = (
            self._n
            * self._fc
            / self._eps_0
            * (1 - (eps[(eps <= 0) & (eps >= self._eps_0)] / self._eps_0))
            ** (self._n - 1)
        )
        # Elsewhere tangent is zero
        tangent[eps < self._eps_0] = 0.0
        tangent[eps > 0] = 0.0
        return tangent

    def __marin__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[float], t.List[float]]:
        """Returns coefficients and strain limits for Marin integration in a
        simply formatted way.

        Arguments:
            strain (float, float): Tuple defining the strain profile: eps =
                strain[0] + strain[1]*y.

        Example:
            [(0, -0.002), (-0.002, -0.003)]
            [(a0, a1, a2), (a0)]
        """
        if self._n != 2:
            # The constitutive law is not writtable as a polynomial,
            # Call the generic distretizing method
            return super().__marin__(strain=strain)

        strains = []
        coeff = []
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # understand in which branch are we
            strain[0] = self.preprocess_strains_with_limits(strain[0])[0]
            if strain[0] > 0:
                # We are in tensile branch
                strains = None
                coeff.append((0.0,))
            elif strain[0] > self._eps_0:
                # We are in the parabolic branch
                strains = None
                a0 = (
                    2
                    * self._fc
                    * strain[0]
                    / self._eps_0
                    * (1 - 0.5 * (strain[0] / self._eps_0))
                )
                a1 = (
                    2
                    * self._fc
                    / self._eps_0
                    * strain[1]
                    * (1 - strain[0] / self._eps_0)
                )
                coeff.append((a0, a1, 0.0))
            elif strain[0] >= self._eps_u:
                # We are in the constant branch
                strains = None
                coeff.append((self._fc,))
            else:
                # We are in a branch of non-resisting concrete
                # Too much compression
                strains = None
                coeff.append((0.0,))
        else:
            # Parabolic part
            strains.append((self._eps_0, 0))
            a0 = (
                2
                * self._fc
                * strain[0]
                / self._eps_0
                * (1 - 0.5 * (strain[0] / self._eps_0))
            )
            a1 = (
                2
                * self._fc
                / self._eps_0
                * strain[1]
                * (1 - strain[0] / self._eps_0)
            )
            a2 = -self._fc * strain[1] ** 2 / self._eps_0**2
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


class BilinearCompression(ConstitutiveLaw):
    """Class for Bilinear Elastic-PerfectlyPlastic Constitutive Law for
    Concrete (only compression behavior).
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_c: float,
        eps_cu: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a BilinearCompression Material.

        Arguments:
            fc (float): Compressive strength (negative number).
            eps_c (float): Strain at compressive strength (pure number).

        Keyword Arguments:
            eps_cu (float): Ultimate strain (pure number).
            name (str): A descriptive name for the constitutive law.
        """
        name = name if name is not None else 'BilinearCompressionLaw'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._eps_c = -abs(eps_c)
        self._eps_cu = -abs(eps_cu)
        self._E = self._fc / self._eps_c

    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Return the stress given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        sig = self._E * eps
        sig[sig < self._fc] = self._fc
        sig[eps > 0] = 0
        sig[eps < self._eps_cu] = 0
        return sig

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent for given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        tangent = np.ones_like(eps) * self._E
        tangent[eps < self._eps_c] = 0.0

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
            # understand in which branch we are
            strain[0] = self.preprocess_strains_with_limits(strain[0])[0]
            if strain[0] > 0:
                # We are in tensile branch
                strains = None
                coeff.append((0.0,))
            elif strain[0] > self._eps_0:
                # We are in the linear branch
                strains = None
                a0 = self._E * strain[0]
                a1 = self._E * strain[1]
                coeff.append((a0, a1))
            elif strain[0] >= self._eps_cu:
                # We are in the constant branch
                strains = None
                coeff.append((self._fc,))
            else:
                # We are in a branch of non-resisting concrete
                # Too much compression
                strains = None
                coeff.append((0.0,))
        else:
            # linear part
            strains.append((self._eps_c, 0))
            a0 = self._E * strain[0]
            a1 = self._E * strain[1]
            coeff.append((a0, a1))
            # Constant part
            strains.append((self._eps_cu, self._eps_c))
            coeff.append((self._fc,))
        return strains, coeff

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        if yielding:
            return (100, self._eps_c)
        return (100, self._eps_cu)


class Sargin(ConstitutiveLaw):
    """Class for Sargin constitutive law.

    The stresses and strains are assumed negative in compression and positive
    in tension.

    References:
    Sargin, M. (1971), "Stress-strain relationship for concrete and the
    analysis of structural concrete section, Study No. 4,
    Solid Mechanics Division, University of Waterloo, Ontario, Canada
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
            fc (float): The strength of concrete in compression

        Keyword Arguments:
            eps_c1 (float): Peak strain of concrete in compression. Default
                value = -0.0023.
            eps_u (float): Ultimate strain of concrete in compression. Default
                value = -0.0035.
            k (float): Plasticity number. Default value = 2.04.
            name (str): A name for the constitutive law.

        Raises:
            ValueError: If k is less or equal to 0.

        Note:
            If positive values are input for fc, eps_c1 and eps_cu1 are input,
            they will be assumed negative.
        """
        name = name if name is not None else 'SarginLaw'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._eps_c1 = -abs(eps_c1)
        self._eps_cu1 = -abs(eps_cu1)
        self._k = k

    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Return the stress given the strain."""
        eps = np.atleast_1d(np.asarray(eps))
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
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


class Popovics(ConstitutiveLaw):
    """Class for Popovics-Mander constitutive law.

    The stresses and strains are assumed negative in compression and positive
    in tension.

    If the relation Ec = 5000 * sqrt(fc) is used for elastic modulus, the
    constitutive law is identical to the one proposed by Mander et al. (1988).

    References:
    Popovics, S., 1973, “A Numerical Approach to the Complete Stress-Strain
    Curve of Concrete”, Cement and Concrete Research, 3(4), 583-599.

    Mander, J.B., Priestley, M.J.N., Park, R., 1988, "Theoretical Stress-Strain
    Model for Confined Concrete", Journal of Structural Engineering, 114(8),
    1804-1826.
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_c: float = -0.002,
        eps_cu: float = -0.0035,
        Ec: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a Popovics Material.

        Arguments:
            fc (float): the strength of concrete in compression

        Keyword Arguments:
            eps_c (float): Peak strain of concrete in compression. Default
                value = -0.002.
            eps_cu (float): Ultimate strain of concrete in compression. Default
                value = -0.0035.
            E (optional float): Elastic modulus of concrete. If None, the
                equation Ec = 5000 * fc**0.5 proposed by Mander et al. (1988)
                is adopted (fc in MPa). Default value = None.
            name (str): A name for the constitutive law.

        Raises:
            ValueError: If E is less or equal to 0.

        Note:
            If positive values are input for fc, eps_c and eps_cu are input,
            they will be assumed negative.
        """
        name = name if name is not None else 'PopovicsLaw'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._eps_c = -abs(eps_c)
        self._eps_cu = -abs(eps_cu)
        if Ec is None:
            # fc in MPa, relation of Mander et al. (1988)
            Ec = 5000 * abs(fc) ** 0.5
        if Ec <= 0:
            raise ValueError('Elastic modulus must be a positive number.')
        E_sec = self._fc / self._eps_c
        self._n = Ec / (Ec - E_sec)

    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Return the stress given the strain."""
        eps = np.atleast_1d(np.asarray(eps))
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        # Compression branch
        eta = eps / self._eps_c

        sig = self._fc * eta * self._n / (self._n - 1 + eta**self._n)

        # Elsewhere stress is 0.0
        sig[eps < self._eps_cu] = 0.0
        sig[eps > 0] = 0.0

        return sig

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compression branch
        eta = eps / self._eps_c

        tangent = (
            (1 - eta**self._n)
            / (self._n - 1 + eta**self._n) ** 2
            * self._n
            * (self._n - 1)
            * self._fc
            / self._eps_c
        )
        # Elsewhere tangent is zero
        tangent[eps < self._eps_cu] = 0.0
        tangent[eps > 0] = 0.0

        return tangent

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        if yielding:
            return (100, self._eps_c)
        return (100, self._eps_cu)


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
        flag: int = 0,
    ) -> None:
        """Initialize a UserDefined constitutive law.

        Arguments:
            x (ArrayLike): Data for strains.
            y (ArrayLike): Data for stresses. Must be of same length as x.

        Keyword Arguments:
            name (Optional, str): A name for the constitutive law.
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

    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Return the stress given strain."""
        eps = np.atleast_1d(np.asarray(eps))
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        return np.interp(eps, self._x, self._y, left=0, right=0)

    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Return the tangent given strain."""
        eps = np.atleast_1d(np.array(eps))

        # Find the segment index for each x value
        indices = np.searchsorted(self._x, eps) - 1

        # Check that indices are within vlaid range
        indices = np.clip(indices, 0, len(self._slopes) - 1)

        # Get the corresponding slopes
        tangent = self._slopes[indices]

        # Elsewhere tangent is zero
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
            strain[0] = self.preprocess_strains_with_limits(strain[0])[0]
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

    def get_ultimate_strain(self, **kwargs) -> t.Tuple[float, float]:
        """Return the ultimate strain (positive and negative)."""
        del kwargs
        return (self._ultimate_strain_p, self._ultimate_strain_n)

    def set_ultimate_strain(
        self, eps_su=t.Union[float, t.Tuple[float, float]]
    ) -> None:
        """Set ultimate strains for Elastic Material if needed.

        Arguments:
            eps_su (float or (float, float)): Defining ultimate strain if a
                single value is provided the same is adopted for both positive
                and negative strains.
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


CONSTITUTIVE_LAWS: t.Dict[str, ConstitutiveLaw] = {
    'elastic': Elastic,
    'elasticplastic': ElasticPlastic,
    'elasticperfectlyplastic': ElasticPlastic,
    'bilinearcompression': BilinearCompression,
    'parabolarectangle': ParabolaRectangle,
    'popovics': Popovics,
    'sargin': Sargin,
}


def get_constitutive_laws_list() -> t.List[str]:
    """Returns a list with valid keywords for constitutive law factory."""
    return list(CONSTITUTIVE_LAWS.keys())


def create_constitutive_law(
    constitutive_law_name: str, material: Material
) -> ConstitutiveLaw:
    """A factory function to create the constitutive law.

    Arguments:
        constitutive_law_name (str): A string defining a valid constitutive law
            type. The available keys can be get with the method
            `get_constitutive_laws_list`.
        material (Material): The material containing the properties needed for
            the definition of the constitutive law.

    Note:
        For working with this facotry function, the material class
        implementations need to provide special dunder methods (__elastic__,
        __parabolarectangle__, etc.) needed for the specific material that
        return the kwargs needed to create the corresponding constitutive
        law object. If the special dunder method is not found an exception
        will be raised.

        If the consitutive law selected is not available for the specific
        material, an exception will be raised.
    """
    law = None
    const_law = CONSTITUTIVE_LAWS.get(constitutive_law_name.lower())
    if const_law is not None:
        method_name = f'__{constitutive_law_name}__'
        # check if the material object has the special method needed
        if hasattr(material, method_name):
            method = getattr(material, method_name)
            if callable(method):
                # get the kwargs from the special dunder method
                kwargs = method()
                # create the constitutive law
                law = const_law(**kwargs)
        else:
            raise ValueError(
                f'Constitutive law {constitutive_law_name} not available for'
                f' material {material.__class__.__name__}'
            )
    else:
        raise ValueError(f'Unknown constitutive law: {constitutive_law_name}')
    return law
