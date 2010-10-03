#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

os.system('rm -rf simulation.dat')
os.system('rm -rf ./plots/simulation_settings.txt')
os.system('rm -rf ./plots/*.pdf')
os.system('rm -rf ./plots/slotted_disc_plots.tar.bz2')

w=3.0
tf=5.0
os.system('./src/simulate {0} {1} > ./plots/simulation_settings.txt'.format(w, tf))

# record is a python file written by the SlottedDiscs C++ class.  If we change
# how we write the records to file, then this python file will be updated
# automatically and should allow the data to be imported correctly.
from record import record_dt

data = np.fromfile('./simulation.dat', dtype=record_dt)

plt.figure()
plt.plot(data[:]['x'], data[:]['y'], 'r-', label='Disc A')
plt.plot(data[:]['cbx'], data[:]['cby'], 'g-', label='Disc B')
plt.plot(data[:]['sox'], data[:]['soy'], 'b-', label='CM')
for i in range(data.size):
    if i % 20 == 0:
        plt.plot([data[i]['x'], data[i]['cbx']], [data[i]['y'], data[i]['cby']],  'k-')
plt.title('Disc contact points')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.savefig('./plots/contactpoints.pdf')

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
plt.title('Reaction forces')
plt.plot(data[:]['t'], data[:]['fay'], 'r--', label=r'$R_a \cdot cl_2$')
plt.plot(data[:]['t'], data[:]['faz'], 'r-', label=r'$R_a \cdot n_3$')
plt.plot(data[:]['t'], data[:]['fby'], 'g--', label=r'$R_b \cdot cl_2$')
plt.plot(data[:]['t'], data[:]['fbz'], 'g-', label=r'$R_b \cdot n_3$')
plt.plot(data[:]['t'], data[:]['fx'], 'b-.', label=r'$(R_a + R_b) \cdot cl_1$')
plt.legend()
plt.xlabel('t [s]')
plt.legend()
plt.savefig('./plots/forces.pdf')


plt.figure()
plt.title('CM acceleration')
plt.plot(data[:]['t'], data[:]['aso1'], label=r'$a_{so} \cdot cl_1$')
plt.plot(data[:]['t'], data[:]['aso2'], label=r'$a_{so} \cdot cl_2$')
plt.plot(data[:]['t'], data[:]['aso3'], label=r'$a_{so} \cdot n_3$')
plt.legend()
plt.xlabel('t [s]')
plt.savefig('./plots/accelerations.pdf')

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
plt.plot(data[:]['t'], data[:]['w'], 'r-', label='w')
plt.xlabel('time [s]')
plt.ylabel('Angular velocity [rad/s]')
plt.title('Angular velocity about contact line')
plt.legend()
plt.savefig('./plots/independentspeed.pdf')

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
plt.subplot(211)
plt.title('Height of COM and Disc B contact')
plt.plot(data[:]['t'], -data[:]['soz'], 'k-', label='cm height')
plt.subplot(212)
plt.plot(data[:]['t'], data[:]['cbz'], 'r-', label='cb height')
plt.xlabel('t [s]')
plt.ylabel(r'meters')
plt.savefig('./plots/heights.pdf')

plt.figure()
plt.title('System angular momentum about mass center')
plt.plot(data[:]['t'], data[:]['H1'], 'r-', label=r'$H^{sys/so}\cdot cl_1$')
plt.plot(data[:]['t'], data[:]['H2'], 'g-', label=r'$H^{sys/so}\cdot cl_2$')
plt.plot(data[:]['t'], data[:]['H3'], 'b-', label=r'$H^{sys/so}\cdot n_3$')
plt.xlabel('t [s]')
plt.ylabel(r'kg m^2 / s')
plt.legend()
plt.savefig('./plots/angularmomentum.pdf')

plt.figure()
plt.title('System linear momentum')
plt.plot(data[:]['t'], data[:]['p1'], 'r-', label=r'$p^{so}\cdot cl_1$')
plt.plot(data[:]['t'], data[:]['p2'], 'g-', label=r'$p^{so}\cdot cl_2$')
plt.plot(data[:]['t'], data[:]['p3'], 'b-', label=r'$p^{so}\cdot n_3$')
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
