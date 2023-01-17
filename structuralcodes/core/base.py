"""Abstract base classes"""
import abc
import typing as t
import warnings


class Material(abc.ABC):
    """Abstract base class for materials."""

    def __init__(self, density: float, name: t.Optional[str] = None) -> None:
        """Initializes an instance of a new material

        Args:
            density (float): density of the material in kg/m3

        Keyword Args:
            name (Optional[str]): descriptive name of the material
        """
        self._density = abs(density)
        self._name = name if name is not None else "Material"

    def update_attributes(self, updated_attributes: t.Dict) -> None:
        """Function for updating the attributes specified in the input
        dictionary

        Args:
            updated_attributes (dict): the dictionary of parameters to be
                updated (not found parameters are skipped with a warning)
        """
        for key, value in updated_attributes.items():
            if not hasattr(self, '_' + key):
                str_list_keys = ''
                for k in updated_attributes.keys():
                    str_list_keys += k + ', '
                str_warn = (
                    f'WARNING: attribute {key} not found. Ignoring the entry.'
                )
                str_warn += '\nAvailable keys: ' + str_list_keys
                warnings.warn(str_warn)
                continue
            setattr(self, '_' + key, value)

    @property
    def name(self):
        """Returns the name of the material"""
        return self._name

    @property
    def density(self):
        """Returns the density of the material in kg/m3"""
        return self._density
