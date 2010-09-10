#!/usr/bin/env python

import matplotlib
import numpy as np
from numpy import sin, cos, abs
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from plotcpp import cpp_record_dt
import os

w3 = 1.0
tf = 2.0
os.system('rm simulatecpp.dat')
os.system('./sim {0} {1} > simulation_settings.txt'.format(w3, tf))
data = np.fromfile('./simulatecpp.dat', dtype=cpp_record_dt)

ra = .1
rb = .1
l = np.sqrt(2.0) * ra
alpha = np.pi/2.0

q2 = np.linspace(-np.pi/2.0, np.pi/2.0, 100)
q3 = np.linspace(-np.pi, np.pi, 200)
Q2, Q3 = np.meshgrid(q2, q3)
hc = rb*(cos(Q2)**2*cos(Q3)**2+(sin(alpha)*sin(Q2)-cos(alpha)*sin(Q3)*cos(Q2))**2)/(1-(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2))**2)**0.5 - l*cos(Q2)*cos(Q3) - ra*cos(Q2)**2/abs(cos(Q2))

plt.figure()
CS = plt.contour(Q2, Q3, hc, levels=[0])
plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel("Lean (q2), radians")
plt.ylabel("Spin (q3), radians")
plt.axis([-np.pi/2.0, np.pi/2.0, -np.pi, np.pi])

plt.hold(True)
plt.plot(data[:]['q2'], data[:]['q3'], 'r.')


plt.savefig('holonomic_constraint.pdf')
