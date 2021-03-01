
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

filename3 = './result/d60.csv'
filename2 = './result/d50.csv'
filename1 = './result/d25.csv'
filename0 = './result/nospray.csv'
# filename1 = './LPT1.csv'
fonts1 = 25
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
CO0 = data0[-2]

data1 = np.loadtxt(filename1, delimiter=',', skiprows=1)
data1 = np.transpose(data1)
grid1 = data1[0]
T1 = data1[1]
U1 = data1[2]
fuel1 = data1[3]
CO1 = data1[-2]

data2 = np.loadtxt(filename2, delimiter=',', skiprows=1)
data2 = np.transpose(data2)
grid2 = data2[0]
T2 = data2[1]
U2 = data2[2]
fuel2 = data2[3]
CO2 = data2[-2]

data3 = np.loadtxt(filename3, delimiter=',', skiprows=1)
data3 = np.transpose(data3)
grid3 = data3[0]
T3 = data3[1]
U3 = data3[2]
fuel3 = data3[3]
CO3 = data3[-2]

#temperature plot:
# plt.scatter(grid0, T0, label='$no spray$',c ='k', s=20)
# plt.scatter(grid1, T1, label='$d = 25$ $\mu m$',c ='r',s=10)
# plt.scatter(grid2, T2, label='$d = 50$ $\mu m$',c ='b',s=10)
# plt.scatter(grid3, T3, label='$d = 60$ $\mu m$',c ='g',s=10)

plt.plot(grid0, T0, label='$Y_{ethanol}$ $no spray$',c ='k')
plt.plot(grid1, T1, label='$d = 25$ $\mu m$',c ='r')
plt.plot(grid2, T2, label='$d = 50$ $\mu m$',c ='b')
plt.plot(grid3, T3, label='$d = 60$ $\mu m$',c ='g')

#Yfuel plot:
# plt.scatter(grid0, fuel0, label='$Y_{ethanol}$ $no spray$',c ='k', s=20)
# plt.scatter(grid1, fuel1, label='$d = 25$ $\mu m$',c ='r',s=10)
# plt.scatter(grid2, fuel2, label='$d = 50$ $\mu m$',c ='b',s=10)
# plt.scatter(grid3, fuel3, label='$d = 60$ $\mu m$',c ='g',s=10)

# plt.plot(grid0, fuel0, label='$Y_{ethanol}$ $no spray$',c ='k')
# plt.plot(grid1, fuel1, label='$d = 25$ $\mu m$',c ='r')
# plt.plot(grid2, fuel2, label='$d = 50$ $\mu m$',c ='b')
# plt.plot(grid3, fuel3, label='$d = 60$ $\mu m$',c ='g')

# #Yco plot:
# plt.scatter(grid0, CO0, label='$Y_{CO}$ $no spray$', marker = 'o', c ='k', s=10)
# plt.scatter(grid1, CO1, label='$d = 25$ $\mu m$', marker = 'o',c ='r',s=6)
# plt.scatter(grid2, CO2, label='$d = 50$ $\mu m$', marker = 'o',c ='b',s=6)
# plt.scatter(grid3, CO3, label='$d = 60$ $\mu m$', marker = 'o',c ='g',s=6)

#temperature plot:
plt.xlabel(r'$z [m]$',fontsize=fonts1)
plt.ylabel(r'$T [K]$',fontsize=fonts1)
#Yfuel plot:
# plt.xlabel(r'$z$ $[m]$',fontsize=fonts1)
# plt.ylabel(r'$Y_{ethanol}$ $[-]$',fontsize=fonts1)
#YCO plot:
# plt.xlabel(r'$z$ $[m]$',fontsize=fonts1)
# plt.ylabel(r'$Y$ $[-]$',fontsize=fonts1)

# ax1.set_ylim(0.0, 1.0)
# ax1.set_xlabel(r'$z$ $(m)$',fontsize=fonts1)
# ax1.set_ylabel(r'$d$ $({\mu}m)$',fontsize=fonts1)

# ax2 = ax1.twinx()
# ax2.scatter(time, Tp/Tp[0], label='Droplet Temperature',c ='b',s=5)
# ax2.set_ylabel(r'$Tp/Tp_0$ $(-)$',fontsize=fonts1)
# ax2.set_ylim(298, 330)

# plt.xlim(0.049, 0.053)
# plt.ylim(290, 2100)
plt.tick_params(labelsize=20)

# handles1, label1 = ax1.get_legend_handles_labels()
# handles2, label2 = ax2.get_legend_handles_labels()
plt.legend(loc='upper right', fontsize=15)
# ax2.legend(loc='upper left', fontsize=12)

plt.show()