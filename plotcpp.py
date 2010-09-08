#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os

os.system('rm simulatecpp.dat')
os.system('./sim')

cpp_record_dt = np.dtype([('t', np.float64),
                          ('q1', np.float64),
                          ('q2', np.float64),
                          ('q3', np.float64),
                          ('x', np.float64),
                          ('y', np.float64),
                          ('w3', np.float64),
                          ('cbx', np.float64),
                          ('cby', np.float64),
                          ('cbz', np.float64),
                          ('ke', np.float64),
                          ('pe', np.float64),
                          ('te', np.float64)])

data = np.fromfile('./simulatecpp.dat', dtype=cpp_record_dt)

plt.figure()
plt.plot(data[:]['t'], data[:]['ke'], 'r-', label='Kinetic')
plt.plot(data[:]['t'], data[:]['pe'], 'g-', label='Potential')
plt.plot(data[:]['t'], data[:]['te'], 'b-', label='Total')
plt.xlabel('t [s]')
plt.ylabel(r'Energy [$kg*m/s^2$]')
plt.title('Energy')
plt.legend()

plt.figure()
plt.plot(data[:]['t'], data[:]['q1'], 'r-', label='Yaw')
plt.plot(data[:]['t'], data[:]['q2'], 'g-', label='Lean')
plt.plot(data[:]['t'], data[:]['q3'], 'b-', label='Spin')
plt.xlabel('t [s]')
plt.ylabel('Angular displacement [rad]')
plt.title('Euler 3-2-3 Angles')
plt.legend()

"""
plt.figure()
plt.plot(data[:]['dcax']-data[0]['dsox'], data[:]['dcay']-data[0]['dsoy'], 'r-', label='SO')
plt.plot(data[0]['dcax']-data[0]['dsox'], data[0]['dcay']-data[0]['dsoy'], 'ro')
plt.plot(data[-1]['dcax']-data[0]['dsox'], data[-1]['dcay']-data[0]['dsoy'], 'rx')
plt.plot(data[:]['x'], data[:]['y'], 'g-', label='SO')
plt.plot(data[0]['x'], data[0]['y'], 'go')
plt.plot(data[-1]['x'], data[-1]['y'], 'gx')
plt.plot(data[:]['dcbx']-data[0]['dsox'], data[:]['dcby']-data[0]['dsoy'], 'b-', label='SO')
plt.plot(data[0]['dcbx']-data[0]['dsox'], data[0]['dcby']-data[0]['dsoy'], 'bo')
plt.plot(data[-1]['dcbx']-data[0]['dsox'], data[-1]['dcby']-data[0]['dsoy'], 'bx')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Position of system center')
plt.legend()

plt.figure()
plt.plot(data[:]['t'], data[:]['w1'], 'r-')
plt.xlabel('t [s]')
plt.ylabel('Angular displacement [rad]')
plt.title('Body fixed angular velocity')
"""

plt.figure()
plt.plot(data[:]['x'], data[:]['y'], 'r-')
plt.plot(data[:]['cbx'], data[:]['cby'], 'g-')
plt.title('Contact points')

plt.figure()
plt.plot(data[:]['t'], data[:]['cbz'], 'r-')
plt.title('Constraints')
plt.show()
