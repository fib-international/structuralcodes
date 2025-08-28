"""Example code for creating an ElasticPlasticMaterial."""

from structuralcodes.materials.basic import ElasticPlasticMaterial

material = ElasticPlasticMaterial(
    E=200000, fy=500, density=7850, Eh=0, eps_su=0.03
)
