import matplotlib.pyplot as plt
import os

def plotcontactpoints(data):
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

def plotangularvelocity(data):
    plt.figure()
    plt.plot(data[:]['t'], data[:]['w1'], 'r-', label='w1')
    plt.plot(data[:]['t'], data[:]['w2'], 'g-', label='w2')
    plt.plot(data[:]['t'], data[:]['w3'], 'b-', label='w3')
    plt.xlabel('t [s]')
    plt.ylabel('Angular velocity [rad/s]')
    plt.title('Body fixed angular velocity')
    plt.legend()
    plt.savefig('./plots/angularvelocity.pdf')

def plotforces(data):
    plt.figure()
    plt.title('Ground reaction forces')
    plt.plot(data[:]['t'], data[:]['fay'], 'r--', label=r'$R_a \cdot cl_2$')
    plt.plot(data[:]['t'], data[:]['faz'], 'r-', label=r'$R_a \cdot n_3$')
    plt.plot(data[:]['t'], data[:]['fby'], 'g--', label=r'$R_b \cdot cl_2$')
    plt.plot(data[:]['t'], data[:]['fbz'], 'g-', label=r'$R_b \cdot n_3$')
    plt.plot(data[:]['t'], data[:]['fx'], 'b-.', label=r'$(R_a + R_b) \cdot cl_1$')
    plt.legend()
    plt.xlabel('t [s]')
    plt.savefig('./plots/forces.pdf')

def plotqdots(data):
    plt.figure()
    plt.title('Coordinate time derivatives')
    plt.plot(data[:]['t'], data[:]['q1p'], 'r-', label=r'$\dot{q}_1$')
    plt.plot(data[:]['t'], data[:]['q2p'], 'g-', label=r'$\dot{q}_2$')
    plt.plot(data[:]['t'], data[:]['q3p'], 'b-', label=r'$\dot{q}_3$')
    plt.plot(data[:]['t'], data[:]['q4p'], 'y-', label=r'$\dot{q}_4$')
    plt.plot(data[:]['t'], data[:]['q5p'], 'k-', label=r'$\dot{q}_5$')
    plt.xlabel('t [s]')
    plt.legend()
    plt.savefig('./plots/qdots.pdf')

def plotaccelerations(data):
    plt.figure()
    plt.title('CM acceleration')
    plt.plot(data[:]['t'], data[:]['aso1'], label=r'$a_{so} \cdot cl_1$')
    plt.plot(data[:]['t'], data[:]['aso2'], label=r'$a_{so} \cdot cl_2$')
    plt.plot(data[:]['t'], data[:]['aso3'], label=r'$a_{so} \cdot n_3$')
    plt.legend()
    plt.xlabel('t [s]')
    plt.savefig('./plots/accelerations.pdf')

def plotenergy(data):
    plt.figure()
    plt.plot(data[:]['t'], data[:]['ke'], 'r-', label='Kinetic')
    plt.plot(data[:]['t'], data[:]['pe'], 'g-', label='Potential')
    plt.plot(data[:]['t'], data[:]['te'], 'b-', label='Total')
    plt.xlabel('t [s]')
    plt.ylabel(r'Energy [kg*m/s^2]')
    plt.title('Energy')
    plt.legend()
    plt.savefig('./plots/energy.pdf')

def plotindependentspeed(data):
    plt.figure()
    plt.plot(data[:]['t'], data[:]['w'], 'r-', label='w')
    plt.xlabel('time [s]')
    plt.ylabel('Angular velocity [rad/s]')
    plt.title('Angular velocity about contact line')
    plt.legend()
    plt.savefig('./plots/independentspeed.pdf')

def ploteulerangles(data):
    plt.figure()
    plt.plot(data[:]['t'], data[:]['q1'], 'r-', label='Yaw')
    plt.plot(data[:]['t'], data[:]['q2'], 'g-', label='Lean')
    plt.plot(data[:]['t'], data[:]['q3'], 'b-', label='Spin')
    plt.xlabel('t [s]')
    plt.ylabel('Angular displacement [rad]')
    plt.title('Euler 3-1-2 Angles')
    plt.legend()
    plt.savefig('./plots/euler312angles.pdf')

def plotheights(data):
    plt.figure()
    plt.subplot(211)
    plt.title('Height of COM and Disc B contact')
    plt.plot(data[:]['t'], -data[:]['soz'], 'k-', label='cm height')
    plt.subplot(212)
    plt.plot(data[:]['t'], data[:]['cbz'], 'r-', label='cb height')
    plt.xlabel('t [s]')
    plt.ylabel(r'meters')
    plt.savefig('./plots/heights.pdf')

def plotangularmomentum(data):
    plt.figure()
    plt.title('System angular momentum about mass center')
    plt.plot(data[:]['t'], data[:]['H1'], 'r-', label=r'$H^{sys/so}\cdot cl_1$')
    plt.plot(data[:]['t'], data[:]['H2'], 'g-', label=r'$H^{sys/so}\cdot cl_2$')
    plt.plot(data[:]['t'], data[:]['H3'], 'b-', label=r'$H^{sys/so}\cdot n_3$')
    plt.xlabel('t [s]')
    plt.ylabel(r'kg m^2 / s')
    plt.legend()
    plt.savefig('./plots/angularmomentum.pdf')

def plotlinearmomentum(data):
    plt.figure()
    plt.title('System linear momentum')
    plt.plot(data[:]['t'], data[:]['p1'], 'r-', label=r'$p^{so}\cdot cl_1$')
    plt.plot(data[:]['t'], data[:]['p2'], 'g-', label=r'$p^{so}\cdot cl_2$')
    plt.plot(data[:]['t'], data[:]['p3'], 'b-', label=r'$p^{so}\cdot n_3$')
    plt.xlabel('t [s]')
    plt.ylabel(r'kg m / s')
    plt.legend()
    plt.savefig('./plots/linearmomentum.pdf')


def plotfunctions(plot_dict, data):
    if plot_dict['contactpoints']:
        plotcontactpoints(data)
    if plot_dict['angularvelocity']:
        plotangularvelocity(data)
    if plot_dict['forces']:
        plotforces(data)
    if plot_dict['qdots']:
        plotqdots(data)
    if plot_dict['accelerations']:
        plotaccelerations(data)
    if plot_dict['energy']:
        plotenergy(data)
    if plot_dict['independentspeed']:
        plotindependentspeed(data)
    if plot_dict['eulerangles']:
        ploteulerangles(data)
    if plot_dict['eulerangles']:
        ploteulerangles(data)
    if plot_dict['heights']:
        plotheights(data)
    if plot_dict['angularmomentum']:
        plotangularmomentum(data)
    if plot_dict['linearmomentum']:
        plotlinearmomentum(data)
    os.system('tar cjf ./plots/slotted_disc_plots.tar.bz2 ./plots/*.pdf' +
              ' ./plots/simulation_settings.txt')

