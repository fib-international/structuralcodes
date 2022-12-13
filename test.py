'''Just a temp test file for this sand-box phase'''
import structuralcodes
from structuralcodes.material.concrete import create_concrete

# We can also create a specific class of material if needed
from structuralcodes.material.concrete import ConcreteMC2010

structuralcodes.set_design_code('mc2010')

# Create a concrete material with the active design code
c25 = create_concrete(25)

print('fcm = ', c25.fcm)
print('fctm = ', c25.fctm)
print('fctkmin = ', c25.fctkmin)

# We can override a quantity (useful for existing materials probably)
c25.fctm = 2.0
print('updated material\n', 'fctm = ', c25.fctm)

# If we want we can reset it to default values combuted by MC2010
c25.fck = 25
print('restored material\n', 'fctm = ', c25.fctm)

# We can also create a material with a given code (not the used default one)
c25 = create_concrete(25, design_code='mc2010')
print('fcm = ', c25.fcm)

c30mc = ConcreteMC2010(30)
print('fcm of the new specific instance: ', c30mc.fcm)

# A new existing concrete object
c30ex = ConcreteMC2010(30, existing=True)
print('Existing concrete - fck: ', c30ex.fck, ' - fcm: ', c30ex.fcm)

# # Compute design strength for new material
# c30 = structuralcodes.Concrete(30)
# print('fcd = ', c30.fcd)

# # Now creae an existing concrete of class C25 (fcm = 33)
# c_existing = structuralcodes.Concrete(25)
# c_existing.existing = True
# print('fcd = ', c_existing.fcd)

# Later on we will be able to do stuff like:
# mySec = structuralcodes.sections.RectangularRCSection(b,h,As,Asprime,c25,
# steel)
# mySec.getMrd('positive')
# chi, M = mySec.computeMomentCurvature('positive')

# and the using sections in structural members to compute other stuff (e.g.
# VRd, etc.)
