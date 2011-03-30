#!/usr/bin/env python

from numpy import sin, cos, arcsin, log, sqrt, linspace, pi
import matplotlib.pyplot as plt

theta = linspace(-2*pi/3+.001, 2*pi/3-.001, 1000)
R = 1.0

st2 = sin(theta/2.0)
ct2 = cos(theta/2.0)
sqrt3 = sqrt(3.0)

xa = 2*R*sqrt3/9*(arcsin(2/sqrt3*st2) +
                  arcsin(st2/sqrt3/ct2) +
                  2.*st2*sqrt(3*ct2*ct2-st2*st2))
ya = 8*R*sqrt3/2*st2*st2-2*R*sqrt3/9*log(ct2)

xb = 2*R*sqrt3/9.0*(arcsin(2.0/sqrt3*st2) +
                    arcsin(st2/sqrt3/ct2) -
                    st2/ct2/ct2*sqrt(3*ct2*ct2-st2*st2))
yb = R*sqrt3/3 + 44*R*sqrt3/9.0*st2*st2 - 2*R*sqrt3/9.0*log(ct2)


plt.plot(xa, ya)
plt.plot(xb, yb)

plt.show()

