"""The concrete class for Model Code 2020 Concrete Material."""

import typing as t
import warnings

from structuralcodes.codes import mc2010

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._concrete import Concrete


class ConcreteMC2010(Concrete):
    """Concrete implementation for MC 2010."""

    # computed values
    _fcm: t.Optional[float] = None
    _fctm: t.Optional[float] = None
    _Eci: t.Optional[float] = None
    _fctkmin: t.Optional[float] = None
    _fctkmax: t.Optional[float] = None
    _Gf: t.Optional[float] = None
    _alpha_cc: t.Optional[float] = None
    _eps_c1: t.Optional[float] = None
    _eps_cu1: t.Optional[float] = None
    _k_sargin: t.Optional[float] = None
    _eps_c2: t.Optional[float] = None
    _eps_cu2: t.Optional[float] = None
    _n_parabolic_rectangular: t.Optional[float] = None
    _eps_c3: t.Optional[float] = None
    _eps_cu3: t.Optional[float] = None

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400.0,
        gamma_c: t.Optional[float] = None,
        alpha_cc: t.Optional[float] = None,
        constitutive_law: t.Optional[
            t.Union[
                t.Literal[
                    'elastic',
                    'parabolarectangle',
                    'bilinearcompression',
                    'sargin',
                    'popovics',
                ],
                ConstitutiveLaw,
            ]
        ] = 'parabolarectangle',
        fcm: t.Optional[float] = None,
        fctm: t.Optional[float] = None,
        fctkmin: t.Optional[float] = None,
        fctkmax: t.Optional[float] = None,
        Eci: t.Optional[float] = None,
        Gf: t.Optional[float] = None,
        eps_c1: t.Optional[float] = None,
        eps_cu1: t.Optional[float] = None,
        k_sargin: t.Optional[float] = None,
        eps_c2: t.Optional[float] = None,
        eps_cu2: t.Optional[float] = None,
        n_parabolic_rectangular: t.Optional[float] = None,
        eps_c3: t.Optional[float] = None,
        eps_cu3: t.Optional[float] = None,
        **kwargs,
    ):
        """Initializes a new instance of Concrete for MC 2010.

        Arguments:
            fck (float): Characteristic strength in MPa if concrete is not
                existing.

        Keyword Arguments:
            name (Optional(str)): A descriptive name for concrete.
            density (float): Density of material in kg/m3 (default: 2400).
            gamma_c (Optional(float)): The partial factor for concrete.
            alpha_cc (float, optional): A factor for considering long-term
                effects on the strength, and effects that arise from the way
                the load is applied.
            fcm (float, optional): The mean compressive strength.
            fctm (float, optional): The mean tensile strength.
            fctkmin (float, optional): The minimum tensile strength.
            fctkmax (float, optional): The maximum tensile strength.
            Eci (float, optional): The initial tangent Young's modulus.
            Gf (float, optional): The tensile fracture energy.
            eps_c1 (float, optional): The strain at peak stress for the Sargin
                constitutive law.
            eps_cu1 (float, optional): The ultimate strain for the Sargin
                constitutive law.
            k_sargin (float, optional): The coefficient for the Sargin
                constitutive law.
            eps_c2 (float, optional): The strain at peak stress for the
                parabolic rectangular constitutive law.
            eps_cu2 (float, optional): The ultimate strain for the parabolic
                rectangular constitutive law.
            n_parabolic_rectangular (float, optional): The coefficient for the
                parabolic rectangular constitutive law.
            eps_c3 (float, optional): The strain at peak stress for the
                bilinear constitutive law.
            eps_cu3 (float, optional): The ultimate strain for the bilinear
                constitutive law.

        Raises:
            ValueError: If fcm is lower than fck.
            ValueError: If k_sargin is negative.
            ValueError: If n_parabolic_rectangular is negative.
            ValueError: If the constitutive law name is not available for the
                material.
            ValueError: If the provided constitutive law is not valid for
                concrete.
            Warning: If Eci is lower than 1e4 or larger than 1e5.
            Warning: If fctm is larger than 0.5 * fck.
            Warning: If eps_c1 is larger than 0.1.
            Warning: If eps_cu1 is larger than 0.1.
            Warning: If eps_c2 is larger than 0.1.
            Warning: If eps_cu2 is larger than 0.1.
            Warning: If n_parabolic_rectangular is larger than 5.
            Warning: If eps_c3 is larger than 0.1.
            Warning: If eps_cu3 is larger than 0.1.
        """
        del kwargs
        if name is None:
            name = f'C{round(fck):d}'
        super().__init__(
            fck=fck,
            name=name,
            density=density,
            gamma_c=gamma_c,
        )
        self._alpha_cc = alpha_cc
        self._fcm = abs(fcm) if fcm is not None else None
        self._fctm = abs(fctm) if fctm is not None else None
        self._fctkmin = abs(fctkmin) if fctkmin is not None else None
        self._fctkmax = abs(fctkmax) if fctkmax is not None else None
        self._Eci = abs(Eci) if Eci is not None else None
        self._Gf = abs(Gf) if Gf is not None else None
        self._eps_c1 = abs(eps_c1) if eps_c1 is not None else None
        self._eps_cu1 = abs(eps_cu1) if eps_cu1 is not None else None
        self._k_sargin = k_sargin if k_sargin is not None else None
        self._eps_c2 = abs(eps_c2) if eps_c2 is not None else None
        self._eps_cu2 = abs(eps_cu2) if eps_cu2 is not None else None
        self._n_parabolic_rectangular = (
            n_parabolic_rectangular
            if n_parabolic_rectangular is not None
            else None
        )
        self._eps_c3 = abs(eps_c3) if eps_c3 is not None else None
        self._eps_cu3 = abs(eps_cu3) if eps_cu3 is not None else None

        self.__post_init__()

        # The constitutive law requires valid attributes, so it should be set
        # after validation
        self._constitutive_law = (
            constitutive_law
            if isinstance(constitutive_law, ConstitutiveLaw)
            else create_constitutive_law(
                constitutive_law_name=constitutive_law, material=self
            )
        )
        if 'concrete' not in self._constitutive_law.__materials__:
            raise ValueError(
                'The provided constitutive law is not valid for concrete.'
            )

    def __post_init__(self):
        """Validator for the attributes that are set in the constructor."""
        # fcm
        if self._fcm is not None and self._fcm <= self._fck:
            raise ValueError(
                (
                    'Mean compressive strength cannot be lower than',
                    'characteristic strength.\n',
                    'Current characteristing strength: ',
                    f'fck = {self._fck}.',
                    f'Current value: value = {self._fcm}',
                )
            )

        # Eci
        if self._Eci is not None and (self._Eci < 1e4 or self._Eci > 1e5):
            warnings.warn(
                'A suspect value of Eci has been input.\n'
                f'Please check Eci that should be in MPa ({self._Eci} given).'
            )

        # fctm
        if self._fctm is not None and self._fctm > 0.5 * self._fck:
            warnings.warn(
                'A suspect value of fctm has been input. Please check.'
            )

        # eps_c1
        if self._eps_c1 is not None and abs(self._eps_c1) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_c1 that should be a pure'
                f' number without units. Please check ({self._eps_c1} given).'
            )

        # eps_cu1
        if self._eps_cu1 is not None and abs(self._eps_cu1) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_cu1 that should be a pure'
                f' number without units. Please check ({self._eps_cu1} given).'
            )

        # k_sargin
        if self._k_sargin is not None and self._k_sargin < 0:
            raise ValueError(
                f'k_sargin should be a positive value ({self._k_sargin} given)'
            )

        # eps_c2
        if self._eps_c2 is not None and abs(self._eps_c2) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_c2 that should be a pure'
                f' number without units. Please check ({self._eps_c2} given).'
            )

        # eps_cu2
        if self._eps_cu2 is not None and abs(self._eps_cu2) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_cu2 that should be a pure'
                f' number without units. Please check ({self._eps_cu2} given).'
            )

        # n_parabolic_rectangular
        if (
            self._n_parabolic_rectangular is not None
            and self._n_parabolic_rectangular < 0
        ):
            raise ValueError(
                'n should be a positive value '
                f'({self._n_parabolic_rectangular} given)'
            )
        if (
            self._n_parabolic_rectangular is not None
            and self._n_parabolic_rectangular >= 5
        ):
            warnings.warn(
                'A suspect value is input for n_parabolic_rectangular. Please '
                'check '
                f'({self._n_parabolic_rectangular} given).'
            )

        # eps_c3
        if self._eps_c3 is not None and abs(self._eps_c3) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_c3 that should be a pure'
                f' number without units. Please check ({self._eps_c3} given).'
            )

        # eps_cu3
        if self._eps_cu3 is not None and abs(self._eps_cu3) >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_cu3 that should be a pure'
                f' number without units. Please check ({self._eps_cu3} given).'
            )

    @property
    def fcm(self) -> float:
        """Returns fcm in MPa.

        Returns:
            float: The mean compressive strength in MPa.

        Note:
            The returned value is derived from fck if fcm is not manually
            provided when initializing the object.
        """
        if self._fcm is None:
            return mc2010.fcm(self._fck)
        return self._fcm

    @property
    def Eci(self) -> float:
        """Returns the modulus of elasticity in MPa at the concrete age of 28
        days.

        It is assumed a normal concrete with quartzite aggregates (alfa_e = 1)

        Returns:
            float: The modulus of elasticity in MPa.

        Note:
            The returned value is derived from fcm if Eci is not manually
            provided when initializing the object.
        """
        if self._Eci is None:
            return mc2010.Eci(self.fcm)
        return self._Eci

    @property
    def fctm(self) -> float:
        """Returns fctm in MPa.

        Returns:
            float: The mean tensile strength in MPa.

        Note:
            The returned value is derived from fck if fctm is not manually
            provided when initializing the object.
        """
        if self._fctm is None:
            return mc2010.fctm(self._fck)
        return self._fctm

    @property
    def fctkmin(self) -> float:
        """Returns fctkmin in MPa.

        Returns:
            float: The lower bound tensile strength in MPa.

        Note:
            The returned value is derived from fctm if fctkmin is not manually
            provided when initializing the object.
        """
        if self._fctkmin is None:
            return mc2010.fctkmin(self.fctm)
        return self._fctkmin

    @property
    def fctkmax(self) -> float:
        """Returns fctkmax in MPa.

        Returns:
            float: The upper bound tensile strength in MPa.

        Note:
            The returned value is derived from fctm if fctkmax is not manually
            provided when initializing the object.
        """
        if self._fctkmax is None:
            return mc2010.fctkmax(self.fctm)
        return self._fctkmax

    @property
    def Gf(self) -> float:
        """Fracture energy of concrete.

        Returns:
            float: The fracture energy in N/m.

        Note:
            The returned value is derived from fck if Gf is not manually
            provided when initializing the object.
        """
        if self._Gf is None:
            return mc2010.Gf(self._fck)
        return self._Gf

    @property
    def gamma_c(self) -> float:
        """The partial factor for concrete."""
        return self._gamma_c or 1.5

    def fcd(self) -> float:
        """Return the design compressive strength in MPa.

        Returns:
            float: The design compressive strength of concrete in MPa.
        """
        # This method should perhaps become a property, but is left as a method
        # for now, to be consistent with other concretes.
        return mc2010.fcd(
            self.fck, alpha_cc=self.alpha_cc, gamma_c=self.gamma_c
        )

    @property
    def alpha_cc(self) -> float:
        """The alpha_cc factor."""
        # Here we should implement the interaction with the globally set
        # national annex. For now, we simply return the default value.
        return self._alpha_cc or 1.0

    @property
    def eps_c1(self) -> float:
        """Returns the strain at maximum compressive strength of concrete (fcm)
        for the Sargin constitutive law.

        Returns:
            float: The strain at maximum compressive strength of concrete.

        Note:
            The returned value is derived from fck if eps_c1 is not manually
            provided when initializing the object.
        """
        if self._eps_c1 is None:
            return mc2010.eps_c1(self._fck)
        return self._eps_c1

    @property
    def eps_cu1(self) -> float:
        """Returns the strain at concrete failure of concrete.

        Returns:
            float: The maximum strength at failure of concrete.

        Note:
            The returned value is derived from fck if eps_cu1 is not manually
            provided when initializing the object.
        """
        if self._eps_cu1 is None:
            return mc2010.eps_cu1(self._fck)
        return self._eps_cu1

    @property
    def k_sargin(self) -> float:
        """Returns the coefficient for Sargin constitutive law.

        Returns:
            float: The plastic coefficient for Sargin law.

        Note:
            The returned value is derived from fck if k_sargin is not manually
            provided when initializing the object.
        """
        if self._k_sargin is None:
            return mc2010.k_sargin(self._fck)
        return self._k_sargin

    @property
    def eps_c2(self) -> float:
        """Returns the strain at maximum compressive strength of concrete (fcd)
        for the Parabola-rectangle constitutive law.

        Returns:
            float: The strain at maximum compressive strength of concrete.

        Note:
            The returned value is derived from fck if eps_c2 is not manually
            provided when initializing the object.
        """
        if self._eps_c2 is None:
            return mc2010.eps_c2(self.fck)
        return self._eps_c2

    @property
    def eps_cu2(self) -> float:
        """Returns the strain at concrete failure of concrete for the
        Parabola-rectangle constitutive law.

        Returns:
            float: The maximum strain at failure of concrete.

        Note:
            The returned value is derived from fck if eps_cu2 is not manually
            provided when initializing the object.
        """
        if self._eps_cu2 is None:
            return mc2010.eps_cu2(self.fck)
        return self._eps_cu2

    @property
    def n_parabolic_rectangular(self) -> float:
        """Returns the coefficient for Parabola-rectangle constitutive law.

        Returns:
            float: The exponent for Parabola-recangle law.

        Note:
            The returned value is derived from fck if n is not manually
            provided when initializing the object.
        """
        if self._n_parabolic_rectangular is None:
            return mc2010.n_parabolic_rectangular(self.fck)
        return self._n_parabolic_rectangular

    @property
    def eps_c3(self) -> float:
        """Returns the strain at maximum compressive strength of concrete (fcd)
        for the Bi-linear constitutive law.

        Returns:
            float: The strain at maximum compressive strength of concrete.

        Note:
            The returned value is derived from fck if eps_c3 is not manually
            provided when initializing the object.
        """
        if self._eps_c3 is None:
            return mc2010.eps_c3(self.fck)
        return self._eps_c3

    @property
    def eps_cu3(self) -> float:
        """Returns the strain at concrete failure of concrete for the Bi-linear
        constitutive law.

        Returns:
            float: The maximum strain at failure of concrete.

        Note:
            The returned value is derived from fck if eps_cu3 is not manually
            provided when initializing the object.
        """
        if self._eps_cu3 is None:
            return mc2010.eps_cu3(self.fck)
        return self._eps_cu3

    def __elastic__(self) -> dict:
        """Returns kwargs for creating an elastic constitutive law."""
        return {'E': self.Eci}

    def __bilinearcompression__(self) -> dict:
        """Returns kwargs for Bi-linear constitutive law."""
        return {
            'fc': self.fcd(),
            'eps_c': self.eps_c3,
            'eps_cu': self.eps_cu3,
        }

    def __parabolarectangle__(self) -> dict:
        """Returns kwargs for creating a parabola rectangle const law."""
        return {
            'fc': self.fcd(),
            'eps_0': self.eps_c2,
            'eps_u': self.eps_cu2,
            'n': self.n_parabolic_rectangular,
        }

    def __sargin__(self) -> dict:
        """Returns kwargs for creating a Sargin const law."""
        return {
            'fc': self.fcd(),
            'eps_c1': self.eps_c1,
            'eps_cu1': self.eps_cu1,
            'k': self.k_sargin,
        }

    def __popovics__(self) -> dict:
        """Returns kwargs for creating a Sargin const law."""
        return {
            'fc': self.fcd(),
            'eps_c': self.eps_c1,
            'eps_cu': self.eps_cu1,
            'Ec': self.Eci,
        }
