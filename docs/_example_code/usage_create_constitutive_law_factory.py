"""Example code for creating a constitutive law using the factory."""

from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.constitutive_laws import create_constitutive_law

# Create a concrete object, and a constitutive law
concrete = ConcreteEC2_2004(fck=45, gamma_c=1.5, alpha_cc=0.85)
constitutive_law = create_constitutive_law(
    constitutive_law_name='bilinearcompression', material=concrete
)
