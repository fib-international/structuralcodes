"""The concrete class for ACI 318-25 Reinforcement Material."""

import typing as t

from structuralcodes.codes import aci318_25

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._reinforcement import Reinforcement


class ReinforcementACI318_25(Reinforcement):  # noqa: N801
    """ACI 318-25 reinforcement material.

    Strengths are unreduced (gamma_s=1.0). ACI applies strength reduction
    factors at member capacity level, not material level.
    Supports ASTM A615 grades: 40, 60, 80, 100.
    """

    def __init__(
        self,
        fyk: float,
        Es: float = 200000.0,
        ftk: float = 550.0,
        epsuk: float = 0.05,
        gamma_s: t.Optional[float] = None,
        name: t.Optional[str] = None,
        density: float = 7850,
        constitutive_law: t.Optional[
            t.Union[
                t.Literal[
                    'elastic',
                    'elasticperfectlyplastic',
                    'elasticplastic',
                ],
                ConstitutiveLaw,
            ]
        ] = 'elasticperfectlyplastic',
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[bool] = None,
        **kwargs,
    ):
        """Initializes a new instance of Reinforcement for ACI 318-25.

        Arguments:
            fyk (float): Characteristic yield strength in MPa.
            Es (float): The Young's modulus in MPa (default: 200000.0).
            ftk (float): Characteristic ultimate strength in MPa
                (default: 550.0).
            epsuk (float): The characteristic strain at the ultimate stress
                level (default: 0.05).
            gamma_s (Optional(float)): The partial factor for reinforcement.
                Default value is 1.0 (ACI applies phi at member level).

        Keyword Arguments:
            name (str): A descriptive name for the reinforcement.
            density (float): Density of material in kg/m3 (default: 7850).
            constitutive_law (ConstitutiveLaw | str): A valid ConstitutiveLaw
                object for reinforcement or a string defining a valid
                constitutive law type for reinforcement. (valid options for
                string: 'elastic', 'elasticplastic', or
                'elasticperfectlyplastic').
            initial_strain (Optional[float]): Initial strain of the material.
            initial_stress (Optional[float]): Initial stress of the material.
            strain_compatibility (Optional[bool]): Only relevant if
                initial_strain or initial_stress are different from zero. If
                True, the material deforms with the geometry. If False, the
                stress in the material upon loading is kept constant
                corresponding to the initial strain.

        Raises:
            ValueError: If the constitutive law name is not available for the
                material.
            ValueError: If the provided constitutive law is not valid for
                reinforcement.
        """
        del kwargs
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

    @classmethod
    def from_grade(
        cls, grade: str = '60', epsuk: float = 0.05, **kwargs
    ) -> 'ReinforcementACI318_25':
        """Create a reinforcement instance from an ASTM A615 grade.

        Arguments:
            grade (str): Reinforcement grade designation. Must be one of
                '40', '60', '80', or '100'. Default is '60'.
            epsuk (float): The characteristic strain at the ultimate stress
                level (default: 0.05).
            **kwargs: Additional keyword arguments passed to the constructor.

        Returns:
            ReinforcementACI318_25: A new reinforcement instance.
        """
        props = aci318_25.reinforcement_grade_props(grade)
        return cls(fyk=props['fy'], ftk=props['fu'], epsuk=epsuk, **kwargs)

    @property
    def gamma_s(self) -> float:
        """The partial factor for reinforcement.

        ACI 318-25 applies strength reduction factors (phi) at the member
        capacity level, not at the material level. The default value is 1.0.
        """
        return self._gamma_s or 1.0

    def fyd(self) -> float:
        """The design yield strength."""
        return self.fyk / self.gamma_s

    def ftd(self) -> float:
        """The design ultimate strength."""
        return self.ftk / self.gamma_s

    def epsud(self) -> float:
        """The design ultimate strain."""
        return self.epsuk

    def __elastic__(self) -> dict:
        """Returns kwargs for creating an elastic constitutive law."""
        return {'E': self.Es}

    def __elasticperfectlyplastic__(self) -> dict:
        """Returns kwargs for ElasticPlastic constitutive law with no strain
        hardening.
        """
        return {
            'E': self.Es,
            'fy': self.fyd(),
            'eps_su': self.epsud(),
        }

    def __elasticplastic__(self) -> dict:
        """Returns kwargs for ElasticPlastic constitutive law with strain
        hardening.
        """
        Eh = (self.ftd() - self.fyd()) / (self.epsud() - self.epsyd)
        return {
            'E': self.Es,
            'fy': self.fyd(),
            'Eh': Eh,
            'eps_su': self.epsud(),
        }
