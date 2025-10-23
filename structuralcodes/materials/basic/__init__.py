"""A collection of basic material classes."""

from ._elastic import ElasticMaterial
from ._elasticplastic import ElasticPlasticMaterial
from ._generic import GenericMaterial

__all__ = [
    'ElasticMaterial',
    'ElasticPlasticMaterial',
    'GenericMaterial',
]
