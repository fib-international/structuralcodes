'''
There should be a module managing several units?
We can:
1. Don't manage them and assume that inside the modules everything is consistent (e.g. N for force, mm for length, s for seconds, t for mass)
2. Manage in a simple way with some variables that converts to the base consistent units (pay attention of variable names)
3. Use some package (e.g. Unum, Pint, units, ) to manage units
'''

mm = 1
N = 1
t = 1

m = 1000 * mm
inch = 25.4 * mm
cm = 10 * mm

kN = 1000 * N
daN = 10 * N
kip = 4448.22 * N

mm2 = 1 * mm * mm
cm2 = 100 * mm2
m2 = 1000000 * mm2
inch2 = 1 * inch * 1 * inch

MPa = 1 * N / (1 * mm2)
Pa = 1 * N / (1 * m2)
ksi = 1 * kip / (1 * inch2)

# a = 100 * mm
# b = 1 * m
# c = a + b