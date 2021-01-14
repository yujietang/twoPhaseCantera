
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

filename0 = './FreeFlame.csv'
filename1 = './FreeFlamewithDrag.csv'
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

data2 = np.loadtxt(filename1, delimiter=',', skiprows=1)
data2 = np.transpose(data2)
time_ = data2[0]
xp_ = data2[1]
d_ = data2[2]
d2_ = data2[3]
Tp_ = data2[5]
T_ = data2[-2]
#Tc = 320K:
time2 = [0.0, 0.00505, 0.00903, 0.013, 0.01608, 0.01998, 0.0239, 0.02809, 0.03194, 0.03592, 0.0399, 0.04396]
value2 = [0.95922, 0.90478, 0.85438, 0.80221, 0.76189, 0.71149, 0.65705, 0.60489, 0.55247, 0.50207, 0.4499, 0.39748]
#Tc = 300K:
# time2 = [0.00116, 0.00608, 0.01113, 0.01605, 0.02208, 0.0291, 0.03514, 0.04111, 0.04715, 0.04996]
# value2 = [0.96971, 0.9222, 0.88453, 0.84884, 0.8033, 0.75382, 0.70631, 0.66273, 0.61916, 0.5975]

ax1  = fig.add_subplot(111)

ax1.scatter(time, d2/d2[0], label='$no$ $Drag$',c ='r',s=20)
ax1.scatter(time_, d2_/d2_[0], label='$Yuen$ $&$ $Chen$ $(1976)$',c ='g', marker = '^', s=10)
# ax1.scatter(time2, value2, label='$OpenFOAM$ $solution$', c = 'k', marker = '+', s=50)
ax1.set_xlabel(r'$t$ $(s)$',fontsize=fonts1)
ax1.set_ylabel(r'$d^2/{d_0}^2$ $(-)$',fontsize=fonts1)
ax1.set_ylim(0.0, 1.0)
# ax1.set_xlabel(r'$z$ $(m)$',fontsize=fonts1)
# ax1.set_ylabel(r'$d$ $({\mu}m)$',fontsize=fonts1)

ax2 = ax1.twinx()
ax2.scatter(time, Tp/Tp[0], label='Droplet Temperature',c ='b',s=5)
ax2.set_ylabel(r'$Tp/Tp_0$ $(-)$',fontsize=fonts1)
# ax2.set_ylim(298, 330)

plt.xlim(0.1425, 0.1445)
plt.tick_params(labelsize=10)

# handles1, label1 = ax1.get_legend_handles_labels()
# handles2, label2 = ax2.get_legend_handles_labels()
ax1.legend(loc='lower left', fontsize=15)
ax2.legend(loc='lower right', fontsize=12)

plt.show()