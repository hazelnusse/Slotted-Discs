#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

os.system('rm -rf simulation.dat')
os.system('rm -rf ./plots/simulation_settings.txt')
os.system('rm -rf ./plots/*.pdf')
os.system('rm -rf ./plots/slotted_disc_plots.tar.bz2')

w=-1.4
tf=5.0
os.system('./src/simulate {0} {1} > ./plots/simulation_settings.txt'.format(w, tf))

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
plt.plot(data[:]['t'], data[:]['w'], 'r-', label='u')
plt.xlabel('t [s]')
plt.ylabel('Angular velocity [rad/s]')
plt.title('Angular velocity about contact line')
plt.legend()
plt.savefig('./plots/independentspeed.pdf')


plt.figure()
plt.plot(data[:]['x'], data[:]['y'], 'r-', label='Disc A')
plt.plot(data[:]['cbx'], data[:]['cby'], 'g-', label='Disc B')
plt.title('Disc contact points')
plt.legend()
plt.savefig('./plots/contactpoints.pdf')

plt.figure()
plt.subplot(211)
plt.title('Height of COM and Disc B contact')
plt.plot(data[:]['t'], data[:]['pe']/9.81/4.0, 'k-', label='cm height')
plt.subplot(212)
plt.plot(data[:]['t'], data[:]['cbz'], 'r-', label='cb height')
plt.xlabel('t [s]')
plt.ylabel(r'meters')
plt.savefig('./plots/heights.pdf')

plt.figure()
plt.title('Angular momentum in contact line frame')
plt.plot(data[:]['t'], data[:]['H1'], 'r-', label='H1')
plt.plot(data[:]['t'], data[:]['H2'], 'g-', label='H2')
plt.plot(data[:]['t'], data[:]['H3'], 'b-', label='H3')
plt.xlabel('t [s]')
plt.ylabel(r'kg m^2 / s')
plt.legend()
plt.savefig('./plots/angularmomentum.pdf')

plt.figure()
plt.title('Linear momentum in contact line frame')
plt.plot(data[:]['t'], data[:]['p1'], 'r-', label='p1')
plt.plot(data[:]['t'], data[:]['p2'], 'g-', label='p2')
plt.plot(data[:]['t'], data[:]['p3'], 'b-', label='p3')
plt.xlabel('t [s]')
plt.ylabel(r'kg m / s')
plt.legend()
plt.savefig('./plots/linearmomentum.pdf')

os.system('tar cjf ./plots/slotted_disc_plots.tar.bz2 ./plots/*.pdf' +
          ' ./plots/simulation_settings.txt')

if len(sys.argv) == 2 and sys.argv[1] == '--silent':
    pass
else:
    plt.show()
