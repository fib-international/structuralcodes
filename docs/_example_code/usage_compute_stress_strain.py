"""Example code for computing stresses with a constitutive law."""

import matplotlib.pyplot as plt
import numpy as np

from structuralcodes.materials.concrete import ConcreteEC2_2004

# Create a concrete object
concrete = ConcreteEC2_2004(fck=45, gamma_c=1.5, alpha_cc=0.85)

# Create an array of strain values based on the ultimate compressive strain,
# and compute corresponding stresses
strains = np.linspace(concrete.constitutive_law.get_ultimate_strain()[0], 0)
stresses = concrete.constitutive_law.get_stress(strains)

# Plot stress-strain relation
fig, ax = plt.subplots()

ax.plot(strains, stresses, '-k')
ax.set_xlim(xmax=0)
ax.set_ylim(ymax=0)
ax.grid()
ax.set_xlabel(r'$\epsilon$ [-]')
ax.set_ylabel(r'$\sigma$ [MPa]')

fig.tight_layout()
fig.savefig('parabola_rectangle.png', dpi=300)
