"""Abstract base classes."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import abc
import typing as t

import numpy as np
from numpy.typing import ArrayLike

import structuralcodes.core._section_results as s_res


class Material(abc.ABC):
    """Abstract base class for materials."""

    _constitutive_law = None
    _initial_strain: t.Optional[float] = None
    _initial_stress: t.Optional[float] = None
    _strain_compatibility: t.Optional[bool] = None

    def __init__(
        self,
        density: float,
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[bool] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initializes an instance of a new material.

        Args:
            density (float): Density of the material in kg/m3.

        Keyword Args:
            initial_strain (Optional[float]): Initial strain of the material.
            initial_stress (Optional[float]): Initial stress of the material.
            strain_compatibility (Optional[bool]): Only relevant if
                initial_strain or initial_stress are different from zero. If
                True, the material deforms with the geometry. If False, the
                stress in the material upon loading is kept constant
                corresponding to the initial strain.
            name (Optional[str]): descriptive name of the material

        Raise:
            ValueError: if both initial_strain and initial_stress are provided
        """
        self._density = abs(density)
        if initial_strain is not None and initial_stress is not None:
            raise ValueError(
                'Both initial_strain and initial_stress cannot be provided.'
            )
        self._initial_strain = initial_strain
        self._initial_stress = initial_stress
        self._strain_compatibility = strain_compatibility
        self._name = name if name is not None else 'Material'

    @property
    def constitutive_law(self) -> ConstitutiveLaw:
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

    @property
    def initial_strain(self):
        """Returns the initial strain of the material."""
        return self._initial_strain

    @property
    def initial_stress(self):
        """Returns the initial stress of the material."""
        return self._initial_stress

    @property
    def strain_compatibility(self):
        """Returns the strain compatibility of the material.

        If true (default), the strain compatibility is enforced
        haveing the same strain as in all other materials of the
        section at the same point. If false, the strain compatibility
        is not enforced and the initial strain is applied to the section
        independently.
        """
        return self._strain_compatibility

    def _apply_initial_strain(self):
        """Wraps the current constitutive law to apply initial strain."""
        strain_compatibility = (
            self._strain_compatibility
            if self._strain_compatibility is not None
            else True
        )
        if self._initial_stress is not None:
            # Specified a stress, compute the strain from it
            self._initial_strain_from_stress()
        if self._initial_strain is not None:
            # Lazy import to avoid circular dependency
            from structuralcodes.materials.constitutive_laws import (  # noqa: PLC0415
                InitialStrain,
            )

            if self._initial_stress is None:
                # Compute the stress from the strain
                self._initial_stress = self._constitutive_law.get_stress(
                    self._initial_strain
                )

            self._constitutive_law = InitialStrain(
                self._constitutive_law,
                self._initial_strain,
                strain_compatibility,
            )

    def _initial_strain_from_stress(self):
        """Computes the initial strain from the initial stress.

        This function is called internally so it assumes that the
        initial stress is not None
        """
        # Iteratively compute the initial strain that gives the desired
        # initial stress. Note that the wrapped law can be nonlinear
        tol = 1e-12
        max_iter = 100
        target_stress = self._initial_stress
        strain = 0.0
        stress = self._constitutive_law.get_stress(strain)
        d_stress = target_stress - stress
        num_iter = 0
        while abs(d_stress) > tol and num_iter < max_iter:
            tangent = self._constitutive_law.get_tangent(strain)
            if tangent == 0:
                raise ValueError(
                    'Tangent modulus = 0 during initial strain computation.'
                )
            d_strain = d_stress / tangent
            strain += d_strain
            stress = self._constitutive_law.get_stress(strain)
            d_stress = target_stress - stress
            num_iter += 1

        if abs(d_stress) > tol:
            raise RuntimeError('Failed to converge for given initial stress.')

        self._initial_strain = strain


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
    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Each constitutive law should provide a method to return the
        stress given the strain level.
        """

    @abc.abstractmethod
    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Each constitutive law should provide a method to return the
        tangent at a given strain level.
        """

    @abc.abstractmethod
    def get_ultimate_strain(self) -> t.Tuple[float, float]:
        """Each constitutive law should provide a method to return the
        ultimate strain (negative and positive).
        """

    def preprocess_strains_with_limits(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Preprocess strain arrays setting those strains sufficiently
        near to ultimate strain limits to exactly ultimate strain limit.
        """
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        eps_min, eps_max = self.get_ultimate_strain()

        if np.isscalar(eps):
            if np.isclose(eps, eps_max, atol=1e-6):
                return eps_max
            if np.isclose(eps, eps_min, atol=1e-6):
                return eps_min
            return eps
        idxs = np.isclose(eps, np.zeros_like(eps) + eps_max, atol=1e-6)
        eps[idxs] = eps_max
        idxs = np.isclose(eps, np.zeros_like(eps) + eps_min, atol=1e-6)
        eps[idxs] = eps_min

        return eps

    def _discretize_law(self) -> ConstitutiveLaw:
        """Discretize the law as a piecewise linear function."""

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

        eps_min, eps_max = self.get_ultimate_strain()
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
        from structuralcodes.materials.constitutive_laws import (  # noqa: PLC0415
            UserDefined,
        )

        return UserDefined(eps, sig)

    def __marin__(self, **kwargs):
        """Function for getting the strain limits and coefficients
        for marin integration.

        By default the law is discretized as a piecewise linear
        function. Then marin coefficients are computed based on this
        discretization.
        """
        piecewise_law = self._discretize_law()
        return piecewise_law.__marin__(**kwargs)

    def __marin_tangent__(self, **kwargs):
        """Function for getting the strain limits and coefficients
        for marin integration of tangent modulus.

        By default the law is discretized as a piecewise linear
        function. Then marin coefficients are computed based on this
        discretization.
        """
        piecewise_law = self._discretize_law()
        return piecewise_law.__marin_tangent__(**kwargs)

    def get_secant(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Method to return the secant at a given strain level."""
        # Adjust eps if it is not scalar
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)

        # Calculate secant for scalar eps
        if np.isscalar(eps):
            if eps != 0:
                sig = self.get_stress(eps)
                return sig / eps
            return self.get_tangent(eps)

        # Calculate secant for array eps
        secant = np.zeros_like(eps)
        strain_is_zero = eps == 0
        strain_is_nonzero = eps != 0
        secant[strain_is_zero] = self.get_tangent(eps[strain_is_zero])
        secant[strain_is_nonzero] = (
            self.get_stress(eps[strain_is_nonzero]) / eps[strain_is_nonzero]
        )
        return secant


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
    def _calculate_gross_section_properties(self) -> s_res.SectionProperties:
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
