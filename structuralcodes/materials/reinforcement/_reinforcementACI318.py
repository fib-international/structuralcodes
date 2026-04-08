"""The reinforcement class for ACI 318-19 Reinforcement Material."""

import typing as t

from structuralcodes.codes import aci318

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._reinforcement import Reinforcement


class ReinforcementACI318(Reinforcement):  # noqa: N801
    """Reinforcement implementation for ACI 318-19.

    Usage philosophy (LRFD vs. partial factor method):
        ACI 318 uses Load and Resistance Factor Design (LRFD), which
        differs from the Eurocode / fib Model Code partial factor
        method used elsewhere in this library. Material strengths
        are not reduced at the material level: gamma_s defaults to
        1.0 and should be left at 1.0 for standard ACI 318 design.

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
        ConcreteACI318 using the 'parabolarectangle' law, pair this
        reinforcement with 'elasticplastic' (the default) or
        'elasticperfectlyplastic'.
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
                reinforcement. Default is 1.0 (ACI does not use
                material partial factors).
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

    def fyd(self) -> float:
        """The design yield strength.

        Note:
            ACI 318 does not reduce material strength. Returns
            fy / gamma_s, which with default gamma_s=1.0 gives
            the unreduced yield strength.
        """
        return aci318.fy_design(self.fyk, phi=1.0 / self.gamma_s)

    @property
    def gamma_s(self) -> float:
        """The partial factor for reinforcement.

        Note:
            Default is 1.0 for ACI 318 (no material partial
            factor).
        """
        return self._gamma_s or 1.0

    def ftd(self) -> float:
        """The design ultimate strength."""
        return self.ftk / self.gamma_s

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
