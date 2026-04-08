"""The concrete class for ACI 318-19 Concrete Material."""

import typing as t

from structuralcodes.codes import aci318

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._concrete import Concrete


class ConcreteACI318(Concrete):  # noqa: N801
    """Concrete implementation for ACI 318-19.

    Note:
        ACI 318 uses specified compressive strength (fc') rather than
        characteristic strength (fck). The parameter is named fck for
        compatibility with the base class, but represents fc' in ACI
        notation. An fc property is provided as an alias.

    Usage philosophy (LRFD vs. partial factor method):
        ACI 318 uses Load and Resistance Factor Design (LRFD), which
        differs from the Eurocode / fib Model Code partial factor
        method used elsewhere in this library. Two consequences for
        users of this class:

        1. Material strengths are not reduced at the material level.
           gamma_c defaults to 1.0 and should be left at 1.0 for
           standard ACI 318 design. Reducing fc' via gamma_c is not
           ACI 318 compliant.

        2. Strength reduction factors (phi, per ACI 318-19 Section
           21.2) must be applied by the user at the member capacity
           level. This class does not apply phi for you.

        As with the Eurocode materials, the user is responsible for
        choosing a constitutive law and parameters suited to the
        limit state being analyzed:

        - Serviceability / cracked-stiffness estimates: use the
          'elastic' law with Ec, parameterized with mean (expected)
          properties.
        - Ultimate limit state (flexure, axial, biaxial interaction
          domains): use 'parabolarectangle' (the default) or
          'bilinearcompression', which apply the ACI 318-19 Section
          22.2 stress block parameters (alpha1=0.85, eps_cu=0.003).
          Compute nominal capacities Pn, Mn with this material, then
          multiply by the appropriate phi factor at the member level
          to obtain phi*Pn and phi*Mn.
    """

    _Ec: t.Optional[float] = None
    _fr: t.Optional[float] = None
    _fct: t.Optional[float] = None
    _wc: float = 2320.0
    _lambda_s: float = 1.0

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2320,
        gamma_c: t.Optional[float] = None,
        constitutive_law: t.Optional[
            t.Union[
                t.Literal[
                    'elastic',
                    'parabolarectangle',
                    'bilinearcompression',
                ],
                ConstitutiveLaw,
            ]
        ] = 'parabolarectangle',
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[bool] = None,
        Ec: t.Optional[float] = None,
        fr: t.Optional[float] = None,
        fct: t.Optional[float] = None,
        wc: float = 2320.0,
        lambda_s: float = 1.0,
        **kwargs,
    ) -> None:
        """Initializes a new instance of Concrete for ACI 318-19.

        Arguments:
            fck (float): Specified compressive strength (fc') in MPa.

        Keyword Arguments:
            name (str): A descriptive name for concrete.
            density (float): Density of material in kg/m3
                (default: 2320).
            gamma_c (float, optional): Partial factor for concrete.
                Default is 1.0 (ACI does not use material partial
                factors).
            constitutive_law (ConstitutiveLaw | str): A valid
                ConstitutiveLaw object or string. Valid options:
                'elastic', 'parabolarectangle',
                'bilinearcompression'.
            initial_strain (Optional[float]): Initial strain.
            initial_stress (Optional[float]): Initial stress.
            strain_compatibility (Optional[bool]): If True, the
                material deforms with the geometry.
            Ec (float, optional): The modulus of elasticity in MPa.
            fr (float, optional): The modulus of rupture in MPa.
            fct (float, optional): The splitting tensile strength
                in MPa.
            wc (float): Unit weight of concrete in kg/m3
                (default: 2320).
            lambda_s (float): Lightweight modification factor
                (default: 1.0 for normalweight).

        Raises:
            ValueError: If the constitutive law is not valid for
                concrete.
        """
        del kwargs
        if name is None:
            name = f'C{round(fck):d}'
        super().__init__(
            fck=fck,
            name=name,
            density=density,
            existing=False,
            gamma_c=gamma_c,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            strain_compatibility=strain_compatibility,
        )
        self._Ec = abs(Ec) if Ec is not None else None
        self._fr = abs(fr) if fr is not None else None
        self._fct = abs(fct) if fct is not None else None
        self._wc = wc
        self._lambda_s = lambda_s

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
        self._apply_initial_strain()

    @property
    def fc(self) -> float:
        """Returns the specified compressive strength (fc') in MPa.

        Returns:
            float: The specified compressive strength in MPa.

        Note:
            This is an alias for fck, using ACI notation.
        """
        return self._fck

    @property
    def gamma_c(self) -> float:
        """The partial factor for concrete.

        Returns:
            float: The partial factor (default 1.0 for ACI 318).

        Note:
            ACI 318 does not use material partial factors. The
            default value of 1.0 maintains compatibility with the
            base class interface.
        """
        return self._gamma_c or 1.0

    @property
    def Ec(self) -> float:
        """Returns the modulus of elasticity in MPa.

        Returns:
            float: The modulus of elasticity in MPa.

        Note:
            Derived from fc and wc if not manually provided.
        """
        if self._Ec is None:
            return aci318.Ec(self._fck, self._wc)
        return self._Ec

    @property
    def fr(self) -> float:
        """Returns the modulus of rupture in MPa.

        Returns:
            float: The modulus of rupture in MPa.

        Note:
            Derived from fc and lambda_s if not manually provided.
        """
        if self._fr is None:
            return aci318.fr(self._fck, self._lambda_s)
        return self._fr

    @property
    def beta1(self) -> float:
        """Returns the Whitney stress block depth factor.

        Returns:
            float: The stress block depth factor (dimensionless).
        """
        return aci318.beta1(self._fck)

    @property
    def eps_cu(self) -> float:
        """Returns the ultimate concrete strain.

        Returns:
            float: The ultimate strain (dimensionless).
        """
        return aci318.eps_cu()

    @property
    def alpha1(self) -> float:
        """Returns the stress block intensity factor.

        Returns:
            float: The stress block intensity (dimensionless).
        """
        return aci318.alpha1()

    @property
    def fct(self) -> float:
        """Returns the splitting tensile strength in MPa.

        Returns:
            float: The splitting tensile strength in MPa.

        Note:
            Derived from fc and lambda_s if not manually provided.
        """
        if self._fct is None:
            return aci318.fct(self._fck, self._lambda_s)
        return self._fct

    def fcd(self) -> float:
        """The design compressive strength in MPa.

        Returns:
            float: The design compressive strength in MPa.

        Note:
            Returns alpha1 * fc / gamma_c. With ACI 318 defaults
            (alpha1=0.85, gamma_c=1.0), this gives 0.85*fc'.
        """
        return self.alpha1 * self._fck / self.gamma_c

    def __elastic__(self) -> dict:
        """Returns kwargs for an elastic constitutive law."""
        return {'E': self.Ec}

    def __parabolarectangle__(self) -> dict:
        """Returns kwargs for a parabola-rectangle constitutive law.

        Note:
            Uses ACI 318 parameters: peak strain at 0.002,
            ultimate strain at 0.003, parabolic exponent n=2.
        """
        return {
            'fc': self.fcd(),
            'eps_0': 0.002,
            'eps_u': self.eps_cu,
            'n': 2,
        }

    def __bilinearcompression__(self) -> dict:
        """Returns kwargs for a bilinear compression constitutive law.

        Note:
            Provides an elastic-perfectly-plastic alternative to
            the parabola-rectangle default. This is not the ACI
            318 Whitney rectangular stress block. The library
            integrates the parabola-rectangle curve directly,
            which ACI 318-19 Section 22.2.2.3 permits as an
            alternative to the Section 22.2.2.4 equivalent
            rectangular stress block simplification.
        """
        return {
            'fc': self.fcd(),
            'eps_c': 0.002,
            'eps_cu': self.eps_cu,
        }
