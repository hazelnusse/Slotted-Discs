#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os

w3 = "1.675"
os.system('./simulate {0}'.format(w3))

record_dt = np.dtype([('t', np.float64),
                      ('yaw', np.float64),
                      ('lean', np.float64),
                      ('spin', np.float64),
                      ('x', np.float64),
                      ('y', np.float64),
                      ('w3', np.float64),
                      ('cbx', np.float64),
                      ('cby', np.float64),
                      ('cbz', np.float64),
                      ('daox', np.float64),
                      ('daoy', np.float64),
                      ('daoz', np.float64),
                      ('dbox', np.float64),
                      ('dboy', np.float64),
                      ('dboz', np.float64),
                      ('ke', np.float64),
                      ('pe', np.float64),
                      ('te', np.float64)])

data = np.fromfile('./simulate.dat', dtype=record_dt)

plt.figure()
plt.plot(data[:]['x'], data[:]['y'], 'r-', label='Disc A')
plt.plot(data[:]['cbx'], data[:]['cby'], 'g-', label='Disc B')
plt.plot(data[0]['x'], data[0]['y'], 'ro')
plt.plot(data[0]['cbx'], data[0]['cby'], 'go')
plt.plot(data[-1]['x'], data[-1]['y'], 'rx')
plt.plot(data[-1]['cbx'], data[-1]['cby'], 'gx')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Position of contact points')
plt.legend()

plt.figure()
plt.plot(data[:]['t'], -data[:]['daoz'], 'r-', label='Disc A')
plt.plot(data[:]['t'], -data[:]['dboz'], 'g-', label='Disc B')
plt.xlabel('t [s]')
plt.ylabel('z [m]')
plt.title('Height of disc centers')
plt.legend()

plt.figure()
plt.plot(data[:]['t'], data[:]['ke'], 'r-', label='Kinetic')
plt.plot(data[:]['t'], data[:]['pe'], 'g-', label='Potential')
plt.plot(data[:]['t'], data[:]['te'], 'b-', label='Total')
plt.xlabel('t [s]')
plt.ylabel(r'Energy [$kg*m/s^2$]')
plt.title('Energy')
plt.legend()

plt.figure()
plt.plot(data[:]['t'], data[:]['w3'], 'r-')
plt.xlabel('t [s]')
plt.ylabel(r'Angular rate [$rad/s$]')
plt.title('Spin rate about line connecting mass centers')

plt.figure()
plt.plot(data[:]['t'], -data[:]['cbz'], 'r-', label='Disc B contact point height')
plt.plot(data[:]['t'], data[:]['te'] - data[0]['te'], 'g-', label='Total' +
         'energy deviation')
plt.xlabel('t [s]')
plt.title('Conserved quantities plot')

plt.figure()
plt.plot(data[:]['t'], data[:]['yaw'], 'r-', label='Yaw')
plt.plot(data[:]['t'], data[:]['lean'], 'g-', label='Lean')
plt.plot(data[:]['t'], data[:]['spin'], 'b-', label='Spin')
plt.xlabel('t [s]')
plt.ylabel('Angular displacement [rad]')
plt.title('Euler 3-1-2 Angles')
plt.legend()

plt.show()
