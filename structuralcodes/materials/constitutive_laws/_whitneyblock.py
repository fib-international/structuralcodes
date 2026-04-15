"""Whitney Block constitutive law for equivalent rectangular stress block."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class WhitneyBlock(ConstitutiveLaw):
    """Equivalent rectangular stress block for section integration.

    This constitutive law represents the equivalent rectangular compressive
    stress distribution used in ACI 318 and other codes (CSA A23.3, AS 3600)
    for computing nominal flexural strength.

    It is NOT a physical stress-strain relationship. It is a code-calibrated
    design idealization that produces the same resultant force and moment as
    the actual nonlinear concrete stress distribution at nominal strength.
    The specific parameters (stress intensity, depth factor, ultimate strain)
    are code-dependent and are supplied by the material class via the
    constitutive law factory pattern (e.g., ConcreteACI318_25.__whitneyblock__()).

    For integration purposes, this is modeled as a piecewise-constant
    stress-strain function. In a linear strain profile with eps_cu at the
    extreme compression fiber:
    - Strain at depth a = beta1*c corresponds to eps_cu*(1-beta1)
    - Stress = fc for strains between eps_cu*(1-beta1) and eps_cu
    - Stress = 0 for strains between 0 and eps_cu*(1-beta1)

    This representation allows both the Marin and Fiber integrators to
    consume the Whitney block without any modification to the section
    analysis pipeline.

    Args:
        fc: Stress block intensity, typically alpha1 * f'c (MPa).
            Stored internally as a negative value (compression).
        beta1: Depth factor mapping neutral axis depth c to block depth
            a = beta1*c.
        eps_cu: Ultimate concrete strain (default 0.003).
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        beta1: float,
        eps_cu: float = 0.003,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a WhitneyBlock constitutive law.

        Arguments:
            fc (float): Stress block intensity, typically alpha1 * f'c.
            beta1 (float): Depth factor mapping neutral axis depth c to
                block depth a = beta1*c.

        Keyword Arguments:
            eps_cu (float): Ultimate concrete strain (default 0.003).
            name (str): A descriptive name for the constitutive law.
        """
        name = name if name is not None else 'WhitneyBlock'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._beta1 = beta1
        self._eps_cu = -abs(eps_cu)
        self._eps_transition = self._eps_cu * (1.0 - beta1)

    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the stress given strain.

        Returns self._fc for strains in the active zone
        [eps_cu, eps_transition], and 0.0 elsewhere.
        """
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        eps = self.preprocess_strains_with_limits(eps=eps)

        if np.isscalar(eps):
            if self._eps_cu <= eps <= self._eps_transition:
                return self._fc
            return 0.0

        sig = np.zeros_like(eps, dtype=float)
        active = (eps >= self._eps_cu) & (eps <= self._eps_transition)
        sig[active] = self._fc
        return sig

    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the tangent for given strain.

        Always 0.0 since the stress block is piecewise constant.
        """
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)

        if np.isscalar(eps):
            return 0.0

        return np.zeros_like(eps, dtype=float)

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        return (self._eps_cu, 0.0)

    def __marin__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
        """Returns coefficients and strain limits for Marin integration.

        Arguments:
            strain (float, float): Tuple defining the strain profile:
                eps = strain[0] + strain[1]*y.
        """
        strains = []
        coeff = []
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            strain[0] = self.preprocess_strains_with_limits(strain[0])
            if self._eps_cu <= strain[0] <= self._eps_transition:
                # In the active (constant stress) zone
                strains = None
                coeff.append((self._fc,))
            else:
                # Outside active zone: zero stress
                strains = None
                coeff.append((0.0,))
        else:
            # Varying strain: two zones
            # Constant stress zone
            strains.append((self._eps_cu, self._eps_transition))
            coeff.append((self._fc,))
            # Zero stress zone
            strains.append((self._eps_transition, 0))
            coeff.append((0.0,))
        return strains, coeff

    def __marin_tangent__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
        """Returns coefficients and strain limits for Marin integration of
        tangent in a simply formatted way.

        Arguments:
            strain (float, float): Tuple defining the strain profile:
                eps = strain[0] + strain[1]*y.
        """
        strains = []
        coeff = []
        if strain[1] == 0:
            # Tangent is always zero
            strains = None
            coeff.append((0.0,))
        else:
            # Single zone covering entire range, tangent is zero
            strains.append((self._eps_cu, 0))
            coeff.append((0.0,))
        return strains, coeff
