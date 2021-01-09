
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

filename0 = './test.csv'
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
time = data1[0]
xp = data1[1]
d = data1[2]
d2 = data1[3]
Tp = data1[5]
T = data1[-2]

ax1  = fig.add_subplot(111)

ax1.scatter(time, d, label='$d$',c ='r',s=1)
ax1.set_xlabel(r'$t [s]$',fontsize=fonts1)
ax1.set_ylabel(r'$d [{\mu}m]$',fontsize=fonts1)
ax1.set_ylim(14, 25)
# ax1.set_xlabel(r'$z$ $(m)$',fontsize=fonts1)
# ax1.set_ylabel(r'$d$ $({\mu}m)$',fontsize=fonts1)

ax2 = ax1.twinx()
ax2.scatter(time, Tp, label='Droplet Temperature',c ='b',s=1)
ax2.set_ylabel(r'$Tp$ $(K)$',fontsize=fonts1)
# ax2.set_ylim(298, 330)

plt.xlim(0.0, 0.05)
plt.tick_params(labelsize=10)

# handles1, label1 = ax1.get_legend_handles_labels()
# handles2, label2 = ax2.get_legend_handles_labels()
ax1.legend(loc='upper left', fontsize=12)
ax2.legend(loc='upper right', fontsize=12)

plt.show()