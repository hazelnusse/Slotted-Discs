#!/usr/bin/env python

import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

ra = rb = 1.0
l = 1.5
alpha = np.pi/2.0

q2 = np.linspace(-np.pi, np.pi, 100)
q3 = np.linspace(-np.pi, np.pi, 100)
Q2, Q3 = np.meshgrid(q2, q3)
hc = (rb*(1.0-(np.cos(alpha)*np.sin(Q2)+np.sin(alpha)*np.cos(Q2)*np.cos(Q3))**2.0)**0.5
    - l*np.cos(Q2)*np.cos(Q3) - ra*np.cos(Q2)**2.0/np.abs(np.cos(Q2)))

plt.figure()
CS = plt.contour(Q2, Q3, hc, levels=[0])
plt.clabel(CS, inline=1, fontsize=10)

plt.show()
