"""Example code for creating an ElasticMaterial from another material."""

from structuralcodes.materials.basic import ElasticMaterial
from structuralcodes.materials.concrete import ConcreteEC2_2004

concrete = ConcreteEC2_2004(fck=45, alpha_cc=0.85, gamma_c=1.5)
elastic_concrete = ElasticMaterial.from_material(concrete)
