import structuralcodes

structuralcodes.set_design_code('mc2010')

c25 =  structuralcodes.Concrete(25)
print('fcm = ',c25.fcm)
print('fctkmin = ',c25.fctkmin)

# Later on we will be able to do stuff like:
# mySec = structuralcodes.sections.RectangularRCSection(b,h,As,Asprime,c25,steel)
# mySec.getMrd('positive')
# chi, M = mySec.computeMomentCurvature('positive')

# and the using sections in structural members to compute other stuff (e.g. VRd, etc.)