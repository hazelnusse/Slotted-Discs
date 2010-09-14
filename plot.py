#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os


w3 = 1.0
tf = 20.0
os.system('rm -rf simulation.dat')
os.system('./src/simulate {0} {1} > ./plots/simulation_settings.txt'.format(w3, tf))

# record is a python file written by the SlottedDiscs C++ class.  If we change
# how we write the records to file, then this python file will be updated
# automatically and should allow the data to be imported correctly.
from record import record_dt

data = np.fromfile('./simulation.dat', dtype=record_dt)

plt.figure()
plt.plot(data[:]['t'], data[:]['ke'], 'r-', label='Kinetic')
plt.plot(data[:]['t'], data[:]['pe'], 'g-', label='Potential')
plt.plot(data[:]['t'], data[:]['te'], 'b-', label='Total')
plt.xlabel('t [s]')
plt.ylabel(r'Energy [kg*m/s^2]')
plt.title('Energy')
plt.legend()
plt.savefig('./plots/energy.pdf')


plt.figure()
plt.plot(data[:]['t'], data[:]['q1'], 'r-', label='Yaw')
plt.plot(data[:]['t'], data[:]['q2'], 'g-', label='Lean')
plt.plot(data[:]['t'], data[:]['q3'], 'b-', label='Spin')
plt.xlabel('t [s]')
plt.ylabel('Angular displacement [rad]')
plt.title('Euler 3-1-2 Angles')
plt.legend()
plt.savefig('./plots/euler312angles.pdf')

plt.figure()
plt.plot(data[:]['t'], data[:]['w1'], 'r-', label='w1')
plt.plot(data[:]['t'], data[:]['w2'], 'g-', label='w2')
plt.plot(data[:]['t'], data[:]['w3'], 'b-', label='w3')
plt.xlabel('t [s]')
plt.ylabel('Angular velocity [rad/s]')
plt.title('Body fixed angular velocity')
plt.legend()
plt.savefig('./plots/angularvelocity.pdf')

plt.figure()
plt.plot(data[:]['x'], data[:]['y'], 'r-', label='Disc A')
plt.plot(data[:]['cbx'], data[:]['cby'], 'g-', label='Disc B')
plt.title('Disc contact points')
plt.legend()
plt.savefig('./plots/contactpoints.pdf')

plt.figure()
plt.plot(data[:]['t'], data[:]['hc'], 'r-', label='Holonomic')
plt.plot(data[:]['t'], data[:]['nh1'], 'g-', label='Nonholonomic 1')
plt.plot(data[:]['t'], data[:]['nh2'], 'b-', label='Nonholonomic 2')
plt.plot(data[:]['t'], data[:]['te']-data[0]['te'], 'k-', label='Change in total energy')
plt.title('Constraints')
plt.legend()
plt.savefig('./plots/constraints.pdf')
plt.show()

os.system('tar cjf ./plots/twindiscplots.tar.bz2 ./plots')

