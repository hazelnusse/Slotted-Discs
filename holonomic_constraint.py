#!/usr/bin/env python

import matplotlib
from numpy import sin, cos, abs, sqrt, fromfile, linspace, pi, meshgrid
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from record import record_dt
import os

w3 = 1.0
tf = 2.0
os.system('rm -rf simulation.dat')
os.system('./src/simulate {0} {1} > simulation_settings.txt'.format(w3, tf))
data = fromfile('./simulation.dat', dtype=record_dt)

ra = .1
rb = .1
l = 0.15 #sqrt(2.0) * ra
alpha = pi/2.0
g = 9.81
m = 4.0
k = 0.075

N = 50
q2 = linspace(-pi/2.0, pi/2.0, N)
q3 = linspace(-pi, pi, 2*N)
Q2, Q3 = meshgrid(q2, q3)

# Holonomic constraint
hc = rb*(cos(Q2)**2*cos(Q3)**2+(sin(alpha)*sin(Q2)-cos(alpha)*sin(Q3)*cos(Q2))**2)/(1-(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2))**2)**0.5 - l*cos(Q2)*cos(Q3) - ra*cos(Q2)**2/abs(cos(Q2))

def pow(b, e):
    return b**e

# Equilibria
equilibria = g*m*(k*rb*sin(Q3)*cos(Q2)*(sin(Q2)-cos(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)))/(pow((1-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)),2)),0.5)*pow((pow((l+ra*cos(Q3)-rb*cos(Q2)*cos(Q3)/pow((1-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)),2)),0.5)),2)+pow(rb,2)*pow((sin(Q2)-cos(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2))),2)/(1-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)),2))+pow((ra*sin(Q3)-rb*(sin(Q3)*cos(Q2)-sin(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)))/pow((1-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)),2)),0.5)),2)),0.5))-sin(Q2)*(ra*sin(Q3)*(l+ra*cos(Q3)-rb*cos(Q2)*cos(Q3)/pow((1-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)),2)),0.5))-(k+ra*cos(Q3))*(ra*sin(Q3)-rb*(sin(Q3)*cos(Q2)-sin(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)))/pow((1-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)),2)),0.5)))/pow((pow((l+ra*cos(Q3)-rb*cos(Q2)*cos(Q3)/pow((1-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)),2)),0.5)),2)+pow(rb,2)*pow((sin(Q2)-cos(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2))),2)/(1-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)),2))+pow((ra*sin(Q3)-rb*(sin(Q3)*cos(Q2)-sin(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)))/pow((1-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)),2)),0.5)),2)),0.5));

plt.figure()
CS1 = plt.contour(Q2, Q3, hc, levels=[0])
plt.clabel(CS1, inline=1, fontsize=10)
plt.xlabel("Lean (q2), radians")
plt.ylabel("Spin (q3), radians")
plt.axis([-pi/2.0, pi/2.0, -pi, pi])
plt.title('Holonomic constraint')
plt.savefig('holonomic_constraint.pdf')

plt.show()

