#!/usr/bin/env python
import matplotlib.pyplot as plt
from numpy import sin, cos, sqrt, fromfile, linspace, pi, meshgrid

ra = .1
rb = .1
l = .15 # sqrt(2.0) * ra
alpha = pi/2.0
g = 9.81
ma = 2.0
mb = 2.0
m = ma + mb
k = l * mb / m

plt.figure()
N = 1000
q2 = linspace(-pi/2.0, pi/2.0, N)
q3 = linspace(-pi, pi, 2*N)
Q2, Q3 = meshgrid(q2, q3)
equilibria = k*rb*sin(Q3)*cos(Q2)*(sin(Q2)-cos(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)))/(pow((1.0-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)), 2.0)),0.5)*pow((pow((l+ra*cos(Q3)-rb*cos(Q2)*cos(Q3)/pow((1.0-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)), 2.0)),0.5)), 2.0)+pow(rb, 2.0)*pow((sin(Q2)-cos(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2))), 2.0)/(1.0-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)), 2.0))+pow((ra*sin(Q3)-rb*(sin(Q3)*cos(Q2)-sin(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)))/pow((1.0-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)), 2.0)),0.5)), 2.0)),0.5)) - sin(Q2)*(ra*sin(Q3)*(l+ra*cos(Q3)-rb*cos(Q2)*cos(Q3)/pow((1.0-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)), 2.0)),0.5))-(k+ra*cos(Q3))*(ra*sin(Q3)-rb*(sin(Q3)*cos(Q2)-sin(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)))/pow((1.0-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)), 2.0)),0.5)))/pow((pow((l+ra*cos(Q3)-rb*cos(Q2)*cos(Q3)/pow((1.0-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)), 2.0)),0.5)), 2.0)+pow(rb, 2.0)*pow((sin(Q2)-cos(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2))), 2.0)/(1.0-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)), 2.0))+pow((ra*sin(Q3)-rb*(sin(Q3)*cos(Q2)-sin(alpha)*(cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)))/pow((1.0-pow((cos(alpha)*sin(Q2)+sin(alpha)*sin(Q3)*cos(Q2)), 2.0)),0.5)), 2.0)),0.5);

CS = plt.contour(Q2, Q3, equilibria, levels=[-0.025, 0.0, .025])
plt.xlabel("Lean (q2), radians")
plt.ylabel("Spin (q3), radians")
plt.axis([-pi/2.0, pi/2.0, -pi, pi])
plt.title('Equilibria')
plt.savefig('./plots/equilibria.pdf')
plt.show()
