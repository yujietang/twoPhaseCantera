import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.size'] = 20
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'

fndc = np.loadtxt('NC12H26.csv',delimiter=',').T
fict = np.loadtxt('IC16H36.csv',delimiter=',').T
ftdc = np.loadtxt('C10H18.csv',delimiter=',').T
ftln = np.loadtxt('C7H8.csv',delimiter=',').T

plt.figure(figsize = (9,6))
plt.plot(fndc[0], fndc[1], label='N-dodecane')
plt.plot(fict[0], fict[1], label='Isocetane')
plt.plot(ftdc[0], ftdc[1], label='Transdecalin')
plt.plot(ftln[0], ftln[1], label='Toluene')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Density (kg/m3)')
plt.savefig('density.png')

plt.figure(figsize = (9,6))
plt.plot(fndc[0], fndc[2], label='N-dodecane')
plt.plot(fict[0], fict[2], label='Isocetane')
plt.plot(ftdc[0], ftdc[2], label='Transdecalin')
plt.plot(ftln[0], ftln[2], label='Toluene')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Vapor pressure (Pa)')
plt.savefig('vapor_pressure.png')

