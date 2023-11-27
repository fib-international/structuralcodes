"""Abstract base classes."""
import abc
import typing as t
import warnings

import structuralcodes.sections._section_results as s_res


class Material(abc.ABC):
    """Abstract base class for materials."""

    __materials__: t.Tuple[str] = ()
    _stress_strain = None

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
        return self._stress_strain

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
    def get_stress(self, eps: float) -> float:
        """Each constitutive law should provide a method to return the
        stress given the strain level.
        """

    @abc.abstractmethod
    def get_tangent(self, eps: float) -> float:
        """Each constitutive law should provide a method to return the
        tangent at a given strain level.
        """

    def __marin__(self):
        """Function for getting the strain limits and coefficients
        for marin integration.
        """
        raise NotImplementedError

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
    ) -> s_res.UltimateBendingMomentResult:
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
