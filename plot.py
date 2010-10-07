#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import plotfunctions as pf

os.system('rm -rf simulation.dat')
os.system('rm -rf ./plots/simulation_settings.txt')
os.system('rm -rf ./plots/*.pdf')
os.system('rm -rf ./plots/slotted_disc_plots.tar.bz2')

w = 3.0
tf = 5.0
displayplots = True

# Allow displaying of plots to be supressed from commandline
if len(sys.argv) == 2:
    if sys.argv[1] == '--silent':
        displayplots = False

os.system('./src/simulate {0} {1} > ./plots/simulation_settings.txt'.format(w, tf))

# record is a python file written by the SlottedDiscs C++ class.  If we change
# how we write the records to file, then this python file will be updated
# automatically and should allow the data to be imported correctly.  
# *** Must be imported after simulation is run ***
from record import record_dt

# Get the data from file and put into a custom data type -- examine
# ./simulation.data for details on all the data fields.
data = np.fromfile('./simulation.dat', dtype=record_dt)

# Dictionary to control which plots are made.  All plots are automatically
# saved to the ./plots subdirectory in .pdf format, and they are also displayed
# on screen by default.
plot_dict = {'contactpoints': True,
             'angularvelocity': True,
             'forces': False,
             'qdots' : False,
             'accelerations': False,
             'energy': False,
             'independentspeed': False,
             'eulerangles': False,
             'heights': False,
             'angularmomentum': False,
             'linearmomentum': False}

# Do the actual plotting
pf.plotfunctions(plot_dict, data)

# Display plots onscreen if that is what is desired
if displayplots:
    plt.show()
