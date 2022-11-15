import structuralcodes

structuralcodes.set_design_code('mc2010')

c25 = structuralcodes.Concrete(25)
print('fcm = ', c25.fcm)
print('fctm = ', c25.fctm)
print('fctkmin = ', c25.fctkmin)

# We can override a quantity (useful for existing materials probably)
c25.fctm = 2.0
print('updated material\n', vars(c25))

# If we want we can reset it to default values combuted by MC2010
c25.fck = 25
print('restored material\n', vars(c25))

# Compute design strength for new material
c30 = structuralcodes.Concrete(30)
print('fcd = ', c30.fcd)

# Now creae an existing concrete of class C25 (fcm = 33)
c_existing = structuralcodes.Concrete(25)
c_existing.existing = True
print('fcd = ', c_existing.fcd)

# Later on we will be able to do stuff like:
# mySec = structuralcodes.sections.RectangularRCSection(b,h,As,Asprime,c25,steel)
# mySec.getMrd('positive')
# chi, M = mySec.computeMomentCurvature('positive')

# and the using sections in structural members to compute other stuff (e.g. VRd, etc.)
