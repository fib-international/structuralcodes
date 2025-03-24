"""Constitutive laws for materials."""

import typing as t

from ...core.base import ConstitutiveLaw, Material
from ._bilinearcompression import BilinearCompression
from ._elastic import Elastic
from ._elasticplastic import ElasticPlastic
from ._parabolarectangle import ParabolaRectangle
from ._popovics import Popovics
from ._sargin import Sargin
from ._userdefined import UserDefined

__all__ = [
    'Elastic',
    'ElasticPlastic',
    'ParabolaRectangle',
    'BilinearCompression',
    'Popovics',
    'Sargin',
    'UserDefined',
    'get_constitutive_laws_list',
    'create_constitutive_law',
]

CONSTITUTIVE_LAWS: t.Dict[str, ConstitutiveLaw] = {
    'elastic': Elastic,
    'elasticplastic': ElasticPlastic,
    'elasticperfectlyplastic': ElasticPlastic,
    'bilinearcompression': BilinearCompression,
    'parabolarectangle': ParabolaRectangle,
    'popovics': Popovics,
    'sargin': Sargin,
}


def get_constitutive_laws_list() -> t.List[str]:
    """Returns a list with valid keywords for constitutive law factory."""
    return list(CONSTITUTIVE_LAWS.keys())


def create_constitutive_law(
    constitutive_law_name: str, material: Material
) -> ConstitutiveLaw:
    """A factory function to create the constitutive law.

    Arguments:
        constitutive_law_name (str): A string defining a valid constitutive law
            type. The available keys can be get with the method
            `get_constitutive_laws_list`.
        material (Material): The material containing the properties needed for
            the definition of the constitutive law.

    Note:
        For working with this facotry function, the material class
        implementations need to provide special dunder methods (__elastic__,
        __parabolarectangle__, etc.) needed for the specific material that
        return the kwargs needed to create the corresponding constitutive
        law object. If the special dunder method is not found an exception
        will be raised.

        If the consitutive law selected is not available for the specific
        material, an exception will be raised.

    Raises:
        ValueError: If the constitutive law is not available for the material.
        ValueError: If the constitutive law name is unknown.
    """
    const_law = CONSTITUTIVE_LAWS.get(constitutive_law_name.lower())
    if const_law is not None:
        method_name = f'__{constitutive_law_name}__'
        # check if the material object has the special method needed
        if hasattr(material, method_name) and callable(
            getattr(material, method_name)
        ):
            method = getattr(material, method_name)
            # get the kwargs from the special dunder method
            kwargs = method()
            # create and return the constitutive law
            return const_law(**kwargs)
        raise ValueError(
            f'Constitutive law {constitutive_law_name} not available for'
            f' material {material.__class__.__name__}'
        )
    raise ValueError(f'Unknown constitutive law: {constitutive_law_name}')
