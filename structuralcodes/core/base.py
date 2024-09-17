"""Abstract base classes."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import abc
import typing as t
import warnings

import numpy as np
from numpy.typing import ArrayLike

import structuralcodes.core._section_results as s_res


class Material(abc.ABC):
    """Abstract base class for materials."""

    _constitutive_law = None

    def __init__(self, density: float, name: t.Optional[str] = None) -> None:
        """Initializes an instance of a new material.

        Args:
            density (float): density of the material in kg/m3

        Keyword Args:
            name (Optional[str]): descriptive name of the material
        """
        self._density = abs(density)
        self._name = name if name is not None else 'Material'

    def update_attributes(self, updated_attributes: t.Dict) -> None:
        """Function for updating the attributes specified in the input
        dictionary.

        Args:
            updated_attributes (dict): the dictionary of parameters to be
                updated (not found parameters are skipped with a warning)
        """
        for key, value in updated_attributes.items():
            if not hasattr(self, '_' + key):
                str_list_keys = ', '.join(updated_attributes.keys())
                str_warn = (
                    f"WARNING: attribute '{key}' not found."
                    " Ignoring the entry.\n"
                    f"Used keys in the call: {str_list_keys}"
                )
                warnings.warn(str_warn)
                continue
            setattr(self, '_' + key, value)

    @property
    def constitutive_law(self):
        """Returns the ConstitutiveLaw of the object."""
        return self._constitutive_law

    @property
    def name(self):
        """Returns the name of the material."""
        return self._name

    @property
    def density(self):
        """Returns the density of the material in kg/m3."""
        return self._density


class ConstitutiveLaw(abc.ABC):
    """Abstract base class for constitutive laws."""

    __materials__: t.Tuple[str] = ()
    constitutive_law_counter: t.ClassVar[int] = 0

    def __init__(self, name: t.Optional[str] = None) -> None:
        """Initialize a ConstitutiveLaw object."""
        self.id = self.constitutive_law_counter
        self._name = name if name is not None else f'ConstitutiveLaw_{self.id}'
        self._increase_global_counter()

    @property
    def name(self):
        """Returns the name of the constitutive law."""
        return self._name

    @classmethod
    def _increase_global_counter(cls):
        cls.constitutive_law_counter += 1

    @abc.abstractmethod
    def get_stress(self, eps: ArrayLike) -> ArrayLike:
        """Each constitutive law should provide a method to return the
        stress given the strain level.
        """

    @abc.abstractmethod
    def get_tangent(self, eps: ArrayLike) -> ArrayLike:
        """Each constitutive law should provide a method to return the
        tangent at a given strain level.
        """

    @abc.abstractmethod
    def get_ultimate_strain(self) -> t.Tuple[float, float]:
        """Each constitutive law should provide a method to return the
        ultimate strain (positive and negative).
        """

    def preprocess_strains_with_limits(self, eps: ArrayLike) -> ArrayLike:
        """Preprocess strain arrays setting those strains sufficiently
        near to ultimate strain limits to exactly ultimate strain limit.
        """
        eps = np.atleast_1d(np.asarray(eps))
        eps_max, eps_min = self.get_ultimate_strain()

        idxs = np.isclose(eps, np.zeros_like(eps) + eps_max, atol=1e-6)
        eps[idxs] = eps_max
        idxs = np.isclose(eps, np.zeros_like(eps) + eps_min, atol=1e-6)
        eps[idxs] = eps_min

        return eps

    def __marin__(self, **kwargs):
        """Function for getting the strain limits and coefficients
        for marin integration.

        By default the law is discretized as a piecewise linear
        function. Then marin coefficients are computed based on this
        discretization.
        """

        # Discretize the constitutive law in a "smart way"
        def find_x_lim(x, y):
            # Check if there are non-zero values for x > 0
            if np.any(y[0:] != 0):
                # Find the last non-zero index for x > 0
                non_zero_indices = np.nonzero(y[0:])[0]
                x_lim_index = 0 + non_zero_indices[-1]
                return x[x_lim_index]
            # All values are zero for x > 0
            return None

        eps_max, eps_min = self.get_ultimate_strain()
        eps_max = min(eps_max, 1)
        # Analise positive branch
        eps = np.linspace(0, eps_max, 10000)
        sig = self.get_stress(eps)
        sig[(sig < np.max(sig) * 1e-6)] = 0
        eps_lim = find_x_lim(eps, sig)
        # Now discretize the function in 10 steps for positive part
        eps_pos = (
            np.linspace(0, -eps_min, 1)
            if eps_lim is None
            else np.linspace(0, eps_lim, 10)
        )
        # Analise negative branch
        eps = np.linspace(0, eps_min, 10000)
        sig = -self.get_stress(eps)
        sig[(sig < np.max(sig) * 1e-6)] = 0
        eps_lim = find_x_lim(-eps, sig)
        # Now discretize the function in 10 steps for negative part
        eps_neg = (
            np.linspace(eps_min, 0, 1, endpoint=False)
            if eps_lim is None
            else np.linspace(-eps_lim, 0, 10, endpoint=False)
        )

        eps = np.concatenate((eps_neg, eps_pos))
        sig = self.get_stress(eps)
        from structuralcodes.materials.constitutive_laws import UserDefined

        # Return Marin coefficients for linearized version
        return UserDefined(eps, sig).__marin__(**kwargs)

    def get_secant(self, eps: float) -> float:
        """Method to return the
        secant at a given strain level.
        """
        if eps != 0:
            sig = self.get_stress(eps)
            return sig / eps
        return self.get_tangent(eps)


class Section(abc.ABC):
    """Abstract base class for a cross secion.
    The section is defined by local axes y and z (mapped to x and y cartesian
    plane respectively).
    """

    section_counter: t.ClassVar[int] = 0

    def __init__(self, name: t.Optional[str] = None) -> None:
        """Initialize a Section object."""
        self.id = self.section_counter
        self._name = name if name is not None else f'Section_{self.id}'
        self._increase_global_counter()

    @property
    def name(self):
        """Returns the name of the section."""
        return self._name

    @classmethod
    def _increase_global_counter(cls):
        cls.section_counter += 1


class SectionCalculator(abc.ABC):
    """Abstract class for SectionCalculators
    defining the interface.
    """

    def __init__(self, section: Section) -> None:
        """Initialization of SectionCalculator object, providing
        a Section object.
        """
        self.section = section

    @abc.abstractmethod
    def _calculate_gross_section_properties(self) -> s_res.GrossProperties:
        """Calculates the gross section properties of the section
        This function is private and called when the section is created
        It stores the result into the result object.

        Returns:
        gross_section_properties (GrossSection)
        """

    @abc.abstractmethod
    def calculate_bending_strength(
        self, theta=0, n=0
    ) -> s_res.UltimateBendingMomentResults:
        """Calculates the bending strength for given inclination of n.a.
        and axial load.

        Input:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n (float, default = 0): axial load applied to the section (+: tension,
        -: compression)

        Return:
        ultimate_bending_moment_result (UltimateBendingMomentResult)
        """

    @abc.abstractmethod
    def calculate_moment_curvature(
        self, theta=0, n=0
    ) -> s_res.MomentCurvatureResults:
        """Calculates the moment-curvature relation for given inclination of
        n.a. and axial load.

        Input:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n (float, default = 0): axial load applied to the section
            (+: tension, -: compression)
        chi_incr (float, default = 1e-8): the curvature increment for the
            analysis

        Return:
        moment_curvature_result (MomentCurvatureResults)
        """
