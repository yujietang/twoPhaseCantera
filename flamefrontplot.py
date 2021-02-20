
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

filename1 = './result/Tg320.csv'
filename0 = './result/flamespeed.csv'
# filename1 = './LPT1.csv'
fonts1 = 18
fonts2 = 25
figs = (13,9)

fig = plt.figure(figsize=figs)
plt.rcParams['font.family'] = 'Serif'

from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'

data0 = np.loadtxt(filename0, delimiter=',', skiprows=1)
data0 = np.transpose(data0)
grid0 = data0[0]
T0 = data0[1]
U0 = data0[2]
fuel0 = data0[3]

data1 = np.loadtxt(filename1, delimiter=',', skiprows=1)
data1 = np.transpose(data1)
grid1 = data1[0]
T1 = data1[1]
U1 = data1[2]
fuel1 = data1[3]

plt.scatter(grid0, T0, label='$no spray$',c ='r',s=5)
plt.scatter(grid1, T1, label='$d = 25 micron$',c ='b',s=5)
plt.xlabel(r'$z [m]$',fontsize=fonts1)
plt.ylabel(r'$T [K]$',fontsize=fonts1)
# ax1.set_ylim(0.0, 1.0)
# ax1.set_xlabel(r'$z$ $(m)$',fontsize=fonts1)
# ax1.set_ylabel(r'$d$ $({\mu}m)$',fontsize=fonts1)

# ax2 = ax1.twinx()
# ax2.scatter(time, Tp/Tp[0], label='Droplet Temperature',c ='b',s=5)
# ax2.set_ylabel(r'$Tp/Tp_0$ $(-)$',fontsize=fonts1)
# ax2.set_ylim(298, 330)

plt.xlim(0.045, 0.055)
plt.tick_params(labelsize=10)

# handles1, label1 = ax1.get_legend_handles_labels()
# handles2, label2 = ax2.get_legend_handles_labels()
plt.legend(loc='upper left', fontsize=12)
# ax2.legend(loc='upper left', fontsize=12)

plt.show()