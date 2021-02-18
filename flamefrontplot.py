
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

filename0 = './flamespeed.csv'
# filename1 = './LPT1.csv'
fonts1 = 18
fonts2 = 25
figs = (13,9)

fig = plt.figure(figsize=figs)
plt.rcParams['font.family'] = 'Serif'

from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'

data1 = np.loadtxt(filename0, delimiter=',', skiprows=1)
data1 = np.transpose(data1)
grid = data1[0]
T = data1[1]
U = data1[2]
fuel = data1[3]

plt.scatter(grid, T, label='$Temperature$',c ='r',s=5)
plt.set_xlabel(r'$z [m]$',fontsize=fonts1)
plt.set_ylabel(r'$T [K]$',fontsize=fonts1)
# ax1.set_ylim(0.0, 1.0)
# ax1.set_xlabel(r'$z$ $(m)$',fontsize=fonts1)
# ax1.set_ylabel(r'$d$ $({\mu}m)$',fontsize=fonts1)

# ax2 = ax1.twinx()
# ax2.scatter(time, Tp/Tp[0], label='Droplet Temperature',c ='b',s=5)
# ax2.set_ylabel(r'$Tp/Tp_0$ $(-)$',fontsize=fonts1)
# ax2.set_ylim(298, 330)

# plt.xlim(0.142, 0.145)
plt.tick_params(labelsize=10)

# handles1, label1 = ax1.get_legend_handles_labels()
# handles2, label2 = ax2.get_legend_handles_labels()
plt.legend(loc='lower left', fontsize=12)
# ax2.legend(loc='upper left', fontsize=12)

plt.show()