"""Example code for creating specific concretes."""

from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.constitutive_laws import Sargin

concrete = ConcreteEC2_2004(fck=45)

concrete = ConcreteEC2_2004(fck=45, fcm=50, constitutive_law='sargin')

constitutive_law = Sargin(fc=53, eps_c1=-0.002, eps_cu1=-0.003, k=2)
concrete = ConcreteEC2_2004(fck=45, constitutive_law=constitutive_law)
