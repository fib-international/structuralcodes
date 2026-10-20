"""Example code for creating material objects with the factory functions."""

from structuralcodes import set_design_code
from structuralcodes.codes import ec2_2004
from structuralcodes.materials.concrete import create_concrete
from structuralcodes.materials.reinforcement import create_reinforcement

fck = 45
fyk = 500
Es = 200000
ductility_class = 'C'

set_design_code('ec2_2004')
concrete = create_concrete(fck=45, alpha_cc=0.85, gamma_c=1.5)
reinforcement = create_reinforcement(
    fyk=fyk,
    Es=Es,
    gamma_s=1.15,
    **ec2_2004.reinforcement_duct_props(
        fyk=fyk, ductility_class=ductility_class
    ),
)
