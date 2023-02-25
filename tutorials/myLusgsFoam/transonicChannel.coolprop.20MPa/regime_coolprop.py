#!/usr/bin/env python3

from CoolProp.CoolProp import PropsSI

pTot = 20e6
Ttot = 400
M2i  = 0.675
gas = 'N2'

H = PropsSI('HMASS', 'P', pTot, 'T', Ttot, gas)
s = PropsSI('SMASS', 'P', pTot, 'T', Ttot, gas)

u = 10.
for iter in range(20):
    h = H - 0.5*u**2
    a = PropsSI('SPEED_OF_SOUND', 'HMASS', h, 'SMASS', s, gas)
    Ma = u/a
    u = M2i*a

p2 = PropsSI('P', 'HMASS', h, 'SMASS', s, gas)
print("p2 = ", p2)
