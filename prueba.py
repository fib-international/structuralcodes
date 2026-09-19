from structuralcodes.codes.aashto_2020 import _deflections
from structuralcodes.codes.aashto_2020 import _punching

mcr = _deflections.Mcr(fc_prime=5, b=200 / 25.4, h=500 / 25.4, d=450 / 25.4)
print(mcr)