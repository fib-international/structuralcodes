"""The reinforcement class for ACI 318-19 Reinforcement Material."""

import typing as t

from structuralcodes.codes import aci318_19

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._reinforcement import Reinforcement


class ReinforcementACI318_19(Reinforcement):  # noqa: N801
    """Reinforcement implementation for ACI 318-19.

    Usage philosophy (LRFD vs. partial factor method):
        ACI 318 uses Load and Resistance Factor Design (LRFD), which
        differs from the Eurocode / fib Model Code partial factor
        method used elsewhere in this library. Material strengths
        are not reduced at the material level: gamma_s is fixed to
        1.0 for standard ACI 318 design. The gamma_s constructor
        argument is accepted only for compatibility with the base
        reinforcement API and cannot be set to a value other than 1.0.

        Strength reduction factors (phi, per ACI 318-19 Section
        21.2) must be applied by the user at the member capacity
        level. This class does not apply phi for you. For example,
        when computing the nominal moment capacity Mn of a
        reinforced section, use this reinforcement with its
        'elasticplastic' or 'elasticperfectlyplastic' constitutive
        law and unreduced fy, then multiply the resulting Mn by the
        appropriate phi factor (e.g. 0.90 for tension-controlled
        flexure) to obtain phi*Mn.

        As with the Eurocode materials, the user chooses the
        constitutive law appropriate to the limit state. For
        ultimate limit state section analysis paired with
        ConcreteACI318_19 using the 'parabolarectangle' law, pair this
        reinforcement with 'elasticplastic' (the default) or
        'elasticperfectlyplastic'.

        The base constructor uses SI units to match the rest of
        structuralcodes: MPa for stress and kg/m3 for density. Use
        from_grade() or from_ksi() when starting from US customary
        ACI inputs.
    """

    def __init__(
        self,
        fyk: float,
        Es: float,
        ftk: float,
        epsuk: float,
        gamma_s: t.Optional[float] = None,
        name: t.Optional[str] = None,
        density: float = 7850.0,
        constitutive_law: t.Optional[
            t.Union[
                t.Literal[
                    'elastic',
                    'elasticperfectlyplastic',
                    'elasticplastic',
                ],
                ConstitutiveLaw,
            ]
        ] = 'elasticplastic',
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[bool] = None,
    ):
        """Initializes a new instance of Reinforcement for ACI 318-19.

        Arguments:
            fyk (float): Specified yield strength (fy) in MPa.
            Es (float): The Young's modulus in MPa.
            ftk (float): Specified ultimate strength (fu) in MPa.
            epsuk (float): The strain at ultimate stress level.

        Keyword Arguments:
            gamma_s (Optional(float)): The partial factor for
                reinforcement. Must be None or 1.0 because ACI does
                not use material partial factors. It is accepted only
                for compatibility with the base reinforcement API;
                non-1.0 values raise ValueError.
            name (str): A descriptive name for the reinforcement.
            density (float): Density in kg/m3 (default: 7850).
            constitutive_law (ConstitutiveLaw | str): A valid
                ConstitutiveLaw or string. Valid options: 'elastic',
                'elasticplastic', 'elasticperfectlyplastic'.
            initial_strain (Optional[float]): Initial strain.
            initial_stress (Optional[float]): Initial stress.
            strain_compatibility (Optional[bool]): If True, the
                material deforms with the geometry.

        Raises:
            ValueError: If the constitutive law is not valid for
                reinforcement.
        """
        if name is None:
            name = f'Reinforcement{round(fyk):d}'

        super().__init__(
            fyk=fyk,
            Es=Es,
            name=name,
            density=density,
            ftk=ftk,
            epsuk=epsuk,
            gamma_s=gamma_s,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            strain_compatibility=strain_compatibility,
        )
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
        if 'steel' not in self._constitutive_law.__materials__:
            raise ValueError(
                'The provided constitutive law is not valid for reinforcement.'
            )
        self._apply_initial_strain()

    def __post_init__(self):
        """Validator for the attributes that are set in the constructor."""
        # gamma_s
        if self._gamma_s is not None and self._gamma_s != 1.0:
            raise ValueError(
                'ACI 318-19 does not use material partial factors. '
                'gamma_s must be 1.0.'
            )

        # Es
        if self._Es <= 0:
            raise ValueError('Es must be positive.')

        # ftk
        if self._ftk < self._fyk:
            raise ValueError(
                'Specified ultimate strength cannot be lower than '
                'specified yield strength.'
            )

        # epsuk
        if self._epsuk <= self.epsyd:
            raise ValueError(
                'Ultimate strain must be larger than yield strain.'
            )

    def fyd(self) -> float:
        """The design yield strength.

        Note:
            ACI 318 does not reduce material strength. Strength
            reduction factors are applied at member level.
        """
        return aci318_19.fy_design(self.fyk)

    @property
    def gamma_s(self) -> float:
        """The partial factor for reinforcement.

        Note:
            ACI 318 does not use material partial factors. The
            value is fixed at 1.0; the constructor accepts gamma_s
            only for base class compatibility and rejects non-1.0
            values.
        """
        return self._gamma_s or 1.0

    def ftd(self) -> float:
        """The design ultimate strength."""
        return self.ftk

    def epsud(self) -> float:
        """The design ultimate strain."""
        return self.epsuk

    def __elastic__(self) -> dict:
        """Returns kwargs for an elastic constitutive law."""
        return {'E': self.Es}

    def __elasticperfectlyplastic__(self) -> dict:
        """Returns kwargs for ElasticPlastic law with no hardening."""
        return {
            'E': self.Es,
            'fy': self.fyd(),
            'eps_su': self.epsud(),
        }

    def __elasticplastic__(self) -> dict:
        """Returns kwargs for ElasticPlastic law with hardening."""
        Eh = (self.ftd() - self.fyd()) / (self.epsud() - self.epsyd)
        return {
            'E': self.Es,
            'fy': self.fyd(),
            'Eh': Eh,
            'eps_su': self.epsud(),
        }

    @classmethod
    def from_grade(
        cls,
        grade: t.Union[str, int],
        epsuk: float,
        gamma_s: t.Optional[float] = None,
        name: t.Optional[str] = None,
        density_pcf: float = 490.0,
        constitutive_law: t.Optional[
            t.Union[
                t.Literal[
                    'elastic',
                    'elasticperfectlyplastic',
                    'elasticplastic',
                ],
                ConstitutiveLaw,
            ]
        ] = 'elasticplastic',
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[bool] = None,
    ) -> 'ReinforcementACI318_19':
        """Create ACI reinforcement from an ASTM grade designation.

        Arguments:
            grade (str | int): ASTM reinforcement grade designation,
                such as '40', '60', '80', or '100'.
            epsuk (float): Strain at ultimate stress level.

        Keyword Arguments:
            density_pcf (float): Material density in lb/ft3. Converted to
                kg/m3 for the base material model. Default is 490 lb/ft3.
            gamma_s (Optional(float)): Must be None or 1.0 because ACI
                does not use material partial factors; non-1.0 values
                raise ValueError.
            name (str): A descriptive name for the reinforcement.

        Returns:
            ReinforcementACI318_19: A reinforcement material storing SI
            values internally.

        Note:
            ACI 318-19 Table 20.2.1.3(a) does not provide a single
            grade-level ultimate strain. The caller must provide epsuk
            based on the applicable reinforcement specification and bar size.
        """
        grade_label = str(grade)
        props = aci318_19.reinforcement_grade_props(grade_label)
        if name is None:
            name = f'Grade {grade_label}'
        return cls(
            fyk=props['fy'],
            Es=aci318_19.Es(),
            ftk=props['fu'],
            epsuk=epsuk,
            gamma_s=gamma_s,
            name=name,
            density=aci318_19.pcf_to_kg_per_m3(density_pcf),
            constitutive_law=constitutive_law,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            strain_compatibility=strain_compatibility,
        )

    @classmethod
    def from_ksi(
        cls,
        fy_ksi: float,
        fu_ksi: float,
        epsuk: float,
        Es_ksi: float = 29000.0,
        gamma_s: t.Optional[float] = None,
        name: t.Optional[str] = None,
        density_pcf: float = 490.0,
        constitutive_law: t.Optional[
            t.Union[
                t.Literal[
                    'elastic',
                    'elasticperfectlyplastic',
                    'elasticplastic',
                ],
                ConstitutiveLaw,
            ]
        ] = 'elasticplastic',
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[bool] = None,
    ) -> 'ReinforcementACI318_19':
        """Create ACI reinforcement from US customary stress inputs.

        Arguments:
            fy_ksi (float): Specified yield strength in ksi.
            fu_ksi (float): Specified ultimate strength in ksi.
            epsuk (float): Strain at ultimate stress level.

        Keyword Arguments:
            Es_ksi (float): Modulus of elasticity in ksi.
                Default is 29000 ksi.
            gamma_s (Optional(float)): Must be None or 1.0 because ACI
                does not use material partial factors; non-1.0 values
                raise ValueError.
            density_pcf (float): Material density in lb/ft3. Converted to
                kg/m3 for the base material model. Default is 490 lb/ft3.

        Returns:
            ReinforcementACI318_19: A reinforcement material storing SI
            values internally.
        """
        if name is None:
            name = f'Reinforcement{round(fy_ksi):d}ksi'
        return cls(
            fyk=aci318_19.ksi_to_mpa(fy_ksi),
            Es=aci318_19.ksi_to_mpa(Es_ksi),
            ftk=aci318_19.ksi_to_mpa(fu_ksi),
            epsuk=epsuk,
            gamma_s=gamma_s,
            name=name,
            density=aci318_19.pcf_to_kg_per_m3(density_pcf),
            constitutive_law=constitutive_law,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            strain_compatibility=strain_compatibility,
        )
