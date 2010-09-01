#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os

os.system('./simulate')

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
plt.plot(data[:]['x'], data[:]['y'], 'r-')
plt.plot(data[:]['cbx'], data[:]['cby'], 'g-')
plt.plot(data[0]['x'], data[0]['y'], 'ro')
plt.plot(data[0]['cbx'], data[0]['cby'], 'go')
plt.plot(data[-1]['x'], data[-1]['y'], 'rx')
plt.plot(data[-1]['cbx'], data[-1]['cby'], 'gx')

plt.figure()
plt.plot(data[:]['t'], -data[:]['daoz'], 'r-')
plt.plot(data[:]['t'], -data[:]['dboz'], 'g-')

plt.figure()
plt.plot(data[:]['t'], data[:]['ke'], 'r-')
plt.plot(data[:]['t'], data[:]['pe'], 'g-')
plt.plot(data[:]['t'], data[:]['te'], 'y-')

plt.figure()
plt.plot(data[:]['t'], data[:]['w3'], 'r-')

plt.figure()
plt.plot(data[:]['t'], data[:]['cbz'], 'r-')
plt.plot(data[:]['t'], data[:]['te'] - data[0]['te'], 'g-')

plt.figure()
plt.plot(data[:]['t'], data[:]['yaw'], 'r-')
plt.plot(data[:]['t'], data[:]['lean'], 'g-')
plt.plot(data[:]['t'], data[:]['spin'], 'b-')

plt.show()
