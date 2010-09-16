#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

os.system('rm -rf simulation.dat')
os.system('rm -rf ./plots/simulation_settings.txt')
os.system('rm -rf ./plots/*.pdf')
os.system('rm -rf ./plots/slotted_disc_plots.tar.bz2')

w3=1.0
tf=20.0
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
plt.xlabel('t [s]')
plt.legend()
plt.savefig('./plots/constraints.pdf')

plt.figure()
plt.plot(data[:]['t'], data[:]['pe']/9.81/4.0, 'k-', label='cm height')
plt.xlabel('t [s]')
plt.ylabel(r'meters')
plt.title('Center of mass height')
plt.savefig('./plots/com_height.pdf')

plt.figure()
plt.subplot(311)
plt.title('Omega expressed in alternate frame')
plt.plot(data[:]['t'], data[:]['w_alt1'], 'r-', label='Along contact line')
plt.legend()
plt.subplot(312)
plt.plot(data[:]['t'], data[:]['w_alt2'], 'g-', label='In plane, normal to contact line')
plt.legend()
plt.subplot(313)
plt.plot(data[:]['t'], data[:]['w_alt3'], 'b-', label='Normal to ground plane')
plt.legend()
plt.xlabel('t [s]')
plt.savefig('./plots/omega_alt.pdf')

os.system('tar cjf ./plots/slotted_disc_plots.tar.bz2 ./plots/*.pdf' +
          ' ./plots/simulation_settings.txt')

if len(sys.argv) == 2 and sys.argv[1] == '--silent':
    pass
else:
    plt.show()
