#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os, sys, time
import plotfunctions as pf


os.system('rm -rf simulation.dat')
os.system('rm -rf ./plots/simulation_settings.txt')
os.system('rm -rf ./plots/*.pdf')
os.system('rm -rf ./plots/slotted_disc_plots.tar.bz2')

w = 8.400090
tf = 1.0
l = 0.1
displayplots = True     # Change to False if you only want the pdf files
output_folder = "./plots/"
simdata = output_folder + "simulation.dat"
sim_settings = output_folder + "simulation_settings.txt"
ts = time.localtime()
t = str(ts[0]) + str(ts[1]) + str(ts[2]) + str(ts[3]) + str(ts[4]) + str(ts[5])
tarball_name = "sim" + t + ".tar.bz2"

# Allow displaying of plots to be supressed from commandline
if len(sys.argv) == 2:
    if sys.argv[1] == '--silent':
        displayplots = False

os.system(('./src/sdsim --omega={0} --time={1} --offset={2} ' +
           '--output={3} > {4}').format(w, tf, l, simdata, sim_settings))

# record is a python file written by the SlottedDiscs C++ class.  Each record
# is comprised of output quantities written sequential to file as double precision
# floating point numbers, at each time step.
# *** Must be imported after simulation is run ***
from record import record_dt

# Get the data from file and put into a custom data type -- examine
# ./simulation.data for details on all the data fields.
data = np.fromfile('./plots/simulation.dat', dtype=record_dt)

# Dictionary to control which plots are generated.  Plots are saved to the
# ./plots subdirectory in pdf format
plot_dict = {'contactpoints': True,
             'angularvelocity': True,
             'forces': True,
             'qdots' : True,
             'accelerations': True,
             'energy': True,
             'independentspeed': True,
             'eulerangles': True,
             'heights': True,
             'angularmomentum': True,
             'linearmomentum': True}

# Do the actual plotting
pf.plotfunctions(plot_dict, data, output_folder)

# Create a tarball of the plots and simulation settings
os.system('mv ./record.py ' + output_folder)
os.system('rm ./record.pyc')
os.system('tar cjf ' + output_folder + tarball_name + " " + 
           output_folder + '*.pdf ' + sim_settings + " " + 
           output_folder + "record.py")

# Cleanup things that were stored in the zip file
os.system('rm ' + output_folder + '*.pdf ' + sim_settings + " " +
          simdata + " " + output_folder + "record.py")

# Display plots onscreen if desired
if displayplots:
    plt.show()
