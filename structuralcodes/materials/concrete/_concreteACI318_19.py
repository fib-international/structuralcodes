"""The concrete class for ACI 318-19 Concrete Material."""

import typing as t

from structuralcodes.codes import aci318_19

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._concrete import Concrete


class ConcreteACI318_19(Concrete):  # noqa: N801
    """Concrete implementation for ACI 318-19.

    Note:
        Use fc for ACI specified compressive strength. The fck parameter is
        still accepted for compatibility with the base class and material
        factory, and is treated as ACI fc' without fractile conversion.

        The base constructor uses SI units to match the rest of
        structuralcodes: MPa for stress and kg/m3 for density.
        The default density is 2400 kg/m3, approximately 150 lb/ft3
        expressed as mass density. It is the base material density,
        not the ACI equilibrium density used in Eq. 19.2.2.1. Use
        wc for that ACI density-dependent Ec and lambda calculation,
        or use from_psi() when starting from US customary inputs.

    Usage philosophy (LRFD vs. partial factor method):
        ACI 318 uses Load and Resistance Factor Design (LRFD), which
        differs from the Eurocode / fib Model Code partial factor
        method used elsewhere in this library. Two consequences for
        users of this class:

        1. Material strengths are not reduced at the material level.
           gamma_c is fixed to 1.0 for standard ACI 318 design.

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
          'bilinearcompression', which use the ACI 318-19 Section
          22.2 stress block parameters (alpha1=0.85, eps_cu=0.003).
          Compute nominal capacities Pn, Mn with this material, then
          multiply by the appropriate phi factor at the member level
          to obtain phi*Pn and phi*Mn.
    """

    _Ec: t.Optional[float] = None
    _fr: t.Optional[float] = None
    _eps_c0: t.Optional[float] = None
    _wc: t.Optional[float] = None
    _lambda_s: t.Optional[float] = None

    def __init__(
        self,
        fc: t.Optional[float] = None,
        name: t.Optional[str] = None,
        density: float = 2400.0,
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
        wc: t.Optional[float] = None,
        lambda_s: t.Optional[float] = None,
        eps_c0: t.Optional[float] = None,
        fck: t.Optional[float] = None,
        **kwargs,
    ) -> None:
        """Initializes a new instance of Concrete for ACI 318-19.

        Arguments:
            fc (float, optional): Specified compressive strength (ACI fc') in
                MPa.

        Keyword Arguments:
            name (str): A descriptive name for concrete.
            density (float): Material density in kg/m3
                (default: 2400). The default is approximately
                150 lb/ft3 expressed as mass density for the base
                material model. It is separate from wc.
            gamma_c (float, optional): Partial factor for concrete.
                Must be None or 1.0 because ACI does not use material
                partial factors. It is accepted only for compatibility
                with the base concrete API; non-1.0 values raise
                ValueError.
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
            wc (float, optional): Equilibrium density of concrete in
                kg/m3 used for ACI density-dependent Ec and lambda.
                If omitted, Ec uses the ACI normalweight expression
                ``4700 * sqrt(fc)`` and lambda is 1.0. This is
                separate from the base material density.
            lambda_s (float, optional): Lightweight modification
                factor. If omitted and wc is provided, it is derived
                from ACI 318-19 Table 19.2.4.1(a).
            eps_c0 (float, optional): Peak concrete strain used by
                parabolic and bilinear constitutive laws.
            fck (float, optional): Compatibility alias for ACI fc' in MPa,
                accepted for the base concrete class and material factory.

        Raises:
            ValueError: If the constitutive law is not valid for
                concrete.
        """
        del kwargs
        if fck is None and fc is None:
            raise ValueError('Either fc or fck must be provided.')
        if fck is not None and fc is not None and fck != fc:
            raise ValueError(
                'fc and fck cannot both be provided with different values.'
            )
        if fck is None:
            fck = fc
        assert fck is not None
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
        self._wc = wc
        self._lambda_s = lambda_s
        self._eps_c0 = abs(eps_c0) if eps_c0 is not None else None

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
        self._apply_initial_strain()

    def __post_init__(self):
        """Validator for the attributes that are set in the constructor."""
        # fc
        if self._fck < aci318_19.MINIMUM_STRESS_BLOCK_FC:
            raise ValueError(
                'Specified compressive strength fc must be at least '
                f'{aci318_19.MINIMUM_STRESS_BLOCK_FC} MPa for ACI 318-19.'
            )

        # gamma_c
        if self._gamma_c is not None and self._gamma_c != 1.0:
            raise ValueError(
                'ACI 318-19 does not use material partial factors. '
                'gamma_c must be 1.0.'
            )

        # wc
        if self._wc is not None and (self._wc < 1440 or self._wc > 2560):
            raise ValueError(
                f'wc={self._wc} must be between 1440 and 2560 kg/m3'
            )

        # lambda_s
        if self._lambda_s is not None and (
            self._lambda_s <= 0 or self._lambda_s > 1.0
        ):
            raise ValueError(
                f'lambda_s={self._lambda_s} must be in the range (0, 1]'
            )

        # eps_c0
        if self._eps_c0 is not None and (
            self._eps_c0 <= 0 or self._eps_c0 >= 0.1
        ):
            raise ValueError(
                'eps_c0 should be a pure number without units. '
                f'Current: {self._eps_c0}.'
            )

    @property
    def fc(self) -> float:
        """Returns the specified compressive strength (fc') in MPa.

        Returns:
            float: The specified compressive strength in MPa.

        Note:
            The inherited fck property returns the same stored value for
            base-class compatibility.
        """
        return self._fck

    @property
    def gamma_c(self) -> float:
        """The partial factor for concrete.

        Returns:
            float: The partial factor (default 1.0 for ACI 318).

        Note:
            ACI 318 does not use material partial factors. The
            value is fixed at 1.0; the constructor accepts gamma_c
            only for base class compatibility and rejects non-1.0
            values.
        """
        return self._gamma_c or 1.0

    @property
    def lambda_s(self) -> float:
        """Returns the lightweight modification factor."""
        if self._lambda_s is not None:
            return self._lambda_s
        if self._wc is None:
            return 1.0
        return aci318_19.lambda_factor(self._wc)

    @property
    def Ec(self) -> float:
        """Returns the modulus of elasticity in MPa.

        Returns:
            float: The modulus of elasticity in MPa.

        Note:
            Derived from fc and wc if not manually provided. If wc is
            omitted, the ACI normalweight expression is used.
        """
        if self._Ec is None:
            return aci318_19.Ec(self._fck, self._wc)
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
            return aci318_19.fr(self._fck, self.lambda_s)
        return self._fr

    @property
    def eps_cu(self) -> float:
        """Returns the ultimate concrete strain.

        Returns:
            float: The ultimate strain (dimensionless).
        """
        return aci318_19.eps_cu()

    @property
    def eps_c0(self) -> float:
        """Returns the peak concrete strain used by constitutive laws.

        Returns:
            float: The peak strain (dimensionless).

        Note:
            ACI 318-19 does not prescribe a complete concrete
            stress-strain curve. The default value is the library's
            peak-strain parameter for the selected constitutive laws.
        """
        if self._eps_c0 is None:
            return aci318_19.eps_c0()
        return self._eps_c0

    @property
    def alpha1(self) -> float:
        """Returns the stress block intensity factor.

        Returns:
            float: The stress block intensity (dimensionless).
        """
        return aci318_19.alpha1()

    @property
    def beta1(self) -> float:
        """Returns the Whitney stress block depth factor.

        Returns:
            float: The stress block depth factor (dimensionless).
        """
        return aci318_19.beta1(self._fck)

    def fcd(self) -> float:
        """The design compressive strength in MPa.

        Returns:
            float: The design compressive strength in MPa.

        Note:
            Returns alpha1 * fc. ACI 318-19 strength reduction factors
            are applied at member level, not material level.
        """
        return self.alpha1 * self._fck

    def __elastic__(self) -> dict:
        """Returns kwargs for an elastic constitutive law."""
        return {'E': self.Ec}

    def __parabolarectangle__(self) -> dict:
        """Returns kwargs for a parabola-rectangle constitutive law.

        Note:
            Uses library stress-strain parameters for ACI 318:
            peak strain from eps_c0, ultimate strain at 0.003,
            and parabolic exponent n=2.
        """
        return {
            'fc': self.fcd(),
            'eps_0': self.eps_c0,
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
            'eps_c': self.eps_c0,
            'eps_cu': self.eps_cu,
        }

    @classmethod
    def from_psi(
        cls,
        fc_psi: float,
        name: t.Optional[str] = None,
        density_pcf: float = 150.0,
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
        Ec_psi: t.Optional[float] = None,
        fr_psi: t.Optional[float] = None,
        wc_pcf: t.Optional[float] = None,
        lambda_s: t.Optional[float] = None,
        eps_c0: t.Optional[float] = None,
    ) -> 'ConcreteACI318_19':
        """Create ACI concrete from US customary inputs.

        Arguments:
            fc_psi (float): Specified compressive strength (ACI fc') in psi.

        Keyword Arguments:
            density_pcf (float): Material density in lb/ft3. Converted to
                kg/m3 for the base material model. Default is 150 lb/ft3.
            Ec_psi (float, optional): Manually specified concrete modulus
                of elasticity in psi. Converted to MPa.
            fr_psi (float, optional): Manually specified modulus of rupture
                in psi. Converted to MPa.
            wc_pcf (float, optional): Equilibrium density in lb/ft3 for
                ACI density-dependent Ec and lambda calculations. Converted
                to kg/m3. If omitted, normalweight Ec and lambda are used.
            lambda_s (float, optional): Lightweight modification factor.
            eps_c0 (float, optional): Peak concrete strain used by
                parabolic and bilinear constitutive laws.

        Returns:
            ConcreteACI318_19: A concrete material storing SI values
            internally.
        """
        if name is None:
            name = f'C{round(fc_psi):d}psi'
        return cls(
            fc=aci318_19.psi_to_mpa(fc_psi),
            name=name,
            density=aci318_19.pcf_to_kg_per_m3(density_pcf),
            gamma_c=gamma_c,
            constitutive_law=constitutive_law,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            strain_compatibility=strain_compatibility,
            Ec=(aci318_19.psi_to_mpa(Ec_psi) if Ec_psi is not None else None),
            fr=(aci318_19.psi_to_mpa(fr_psi) if fr_psi is not None else None),
            wc=(
                aci318_19.pcf_to_kg_per_m3(wc_pcf)
                if wc_pcf is not None
                else None
            ),
            lambda_s=lambda_s,
            eps_c0=eps_c0,
        )
