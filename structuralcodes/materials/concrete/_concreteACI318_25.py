"""The concrete class for ACI 318-25 Concrete Material."""

import typing as t

from structuralcodes.codes import aci318_25

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._concrete import Concrete


class ConcreteACI318_25(Concrete):  # noqa: N801
    """ACI 318-25 concrete material.

    Uses LRFD philosophy -- material strengths are unreduced. Safety is applied
    at member capacity level via strength reduction factors phi (Ch. 21).

    The gamma_c property returns 1.0 to satisfy the Concrete base class
    interface. ACI 318 does not use material partial factors. The fcd() method
    returns alpha1 * f'c (= 0.85 * f'c), which is the stress intensity used
    in the Whitney equivalent rectangular stress block, not a gamma-reduced
    design strength in the Eurocode sense.
    """

    _Ec: t.Optional[float] = None
    _fr: t.Optional[float] = None
    _wc: float = 2320.0
    _lambda_s: float = 1.0

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400,
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
        wc: float = 2320.0,
        lambda_s: float = 1.0,
        **kwargs,
    ) -> None:
        """Initializes a new instance of Concrete for ACI 318-25.

        Arguments:
            fck (float): Specified compressive strength f'c in MPa.

        Keyword Arguments:
            name (str): A descriptive name for concrete.
            density (float): Density of material in kg/m3 (default: 2400).
            gamma_c (float, optional): Partial factor for concrete. ACI 318
                does not use material partial factors; defaults to 1.0.
            constitutive_law (ConstitutiveLaw | str): A valid ConstitutiveLaw
                object for concrete or a string defining a valid constitutive
                law type for concrete. (valid options for string: 'elastic',
                'parabolarectangle', 'bilinearcompression').
            initial_strain (Optional[float]): Initial strain of the material.
            initial_stress (Optional[float]): Initial stress of the material.
            strain_compatibility (Optional[bool]): Only relevant if
                initial_strain or initial_stress are different from zero. If
                True, the material deforms with the geometry. If False, the
                stress in the material upon loading is kept constant
                corresponding to the initial strain.
            Ec (float, optional): Modulus of elasticity in MPa. If not
                provided, computed per ACI 318-25 Table 19.2.2.1.
            fr (float, optional): Modulus of rupture in MPa. If not provided,
                computed per ACI 318-25 Eq. 19.2.3.1.
            wc (float): Unit weight of concrete in kg/m3 (default: 2320).
            lambda_s (float): Lightweight modification factor (default: 1.0).

        Raises:
            ValueError: If the constitutive law name is not available for the
                material.
            ValueError: If the provided constitutive law is not valid for
                concrete.
            ValueError: If the constitutive law name is unknown.
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
        self._wc = wc
        self._lambda_s = lambda_s

        # The constitutive law requires valid attributes, so it should be set
        # after storing properties
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
        """Returns f'c in MPa (ACI notation alias for fck).

        Returns:
            float: The specified compressive strength in MPa.
        """
        return self._fck

    @property
    def gamma_c(self) -> float:
        """The partial factor for concrete.

        ACI 318 does not use material partial factors. Returns 1.0 by default
        to satisfy the Concrete base class interface.
        """
        return self._gamma_c or 1.0

    @property
    def alpha1(self) -> float:
        """Stress block intensity factor per ACI 318-25, Section 22.2.2.4.1.

        Returns:
            float: alpha1 = 0.85.
        """
        return aci318_25.alpha1()

    def fcd(self) -> float:
        """Return the design compressive strength in MPa.

        For ACI 318-25 this is alpha1 * f'c / gamma_c. Since gamma_c = 1.0,
        this is effectively alpha1 * f'c = 0.85 * f'c, which is the stress
        intensity used in the Whitney equivalent rectangular stress block.

        Returns:
            float: The design compressive strength of concrete in MPa.
        """
        return self.alpha1 * self.fc / self.gamma_c

    @property
    def Ec(self) -> float:
        """Returns Ec in MPa.

        Returns:
            float: The modulus of elasticity of concrete in MPa.

        Note:
            The returned value is computed per ACI 318-25 Table 19.2.2.1 if
            Ec is not manually provided when initializing the object.
        """
        if self._Ec is None:
            return aci318_25.Ec(self.fc, wc=self._wc)
        return self._Ec

    @property
    def fr(self) -> float:
        """Returns the modulus of rupture in MPa.

        Returns:
            float: The modulus of rupture fr in MPa.

        Note:
            The returned value is computed per ACI 318-25 Eq. 19.2.3.1 if
            fr is not manually provided when initializing the object.
        """
        if self._fr is None:
            return aci318_25.fr(self.fc, lambda_s=self._lambda_s)
        return self._fr

    @property
    def fct(self) -> float:
        """Returns the splitting tensile strength in MPa.

        Returns:
            float: The splitting tensile strength fct in MPa.
        """
        return aci318_25.fct(self.fc, lambda_s=self._lambda_s)

    @property
    def beta1(self) -> float:
        """Whitney stress block depth factor per ACI 318-25, Table 22.2.2.4.3.

        Returns:
            float: Stress block depth factor beta1 (dimensionless).
        """
        return aci318_25.beta1(self.fc)

    @property
    def eps_cu(self) -> float:
        """Maximum usable compressive strain in concrete.

        ACI 318-25, Section 22.2.2.1.

        Returns:
            float: Ultimate concrete strain (dimensionless), equal to 0.003.
        """
        return aci318_25.eps_cu()

    def __elastic__(self) -> dict:
        """Returns kwargs for creating an elastic constitutive law."""
        return {'E': self.Ec}

    def __parabolarectangle__(self) -> dict:
        """Returns kwargs for creating a parabola rectangle const law."""
        return {
            'fc': self.fcd(),
            'eps_0': 0.002,
            'eps_u': self.eps_cu,
            'n': 2,
        }

    def __bilinearcompression__(self) -> dict:
        """Returns kwargs for Bi-linear constitutive law."""
        return {
            'fc': self.fcd(),
            'eps_c': 0.002,
            'eps_cu': self.eps_cu,
        }
