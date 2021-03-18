
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

filename0 = './result/2way.csv'
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
z = data1[0]
ug = data1[1]
rhog = data1[2]
Tg = data1[3]
Yf = data1[4]
phig = data1[5]
Sm = data1[-3]
Sh = data1[-2]
Syf = data1[-1]

ax1  = fig.add_subplot(111)

# ax1.scatter(xp, d2/d2[0], label='$d$',c ='r',s=5)
# ax1.scatter(xp, d, label='$d$',c ='r',s=5)
# ax1.scatter(z, Sm, label='${S_m}$', c='r', s=3)
ax1.plot(z, Tg, label='${S_m}$', c='r')
ax1.set_xlabel(r'$z [m]$',fontsize=fonts1)
# ax1.set_ylabel(r'$d^2/{d_0}^2 (-)$',fontsize=fonts1)
ax1.set_ylabel(r'${S_m}$ $({kg/m3s})$',fontsize=fonts1)

# ax1.set_ylim(-2e-5, 0.0)
# ax1.set_xlabel(r'$z$ $(m)$',fontsize=fonts1)
# ax1.set_ylabel(r'$d$ $({\mu}m)$',fontsize=fonts1)

ax2 = ax1.twinx()
# ax2.scatter(xp, Tp/Tp[0], label='Droplet Temperature',c ='b',s=5)
# ax2.scatter(xp, Tp, label='Droplet Temperature',c ='b',s=5)
# ax2.scatter(z, Sh, label='${S_h}$',c ='b',s=3)
ax2.plot(z, Yf, label='${S_h}$',c ='b')

ax2.set_ylabel(r'${S_h}$ $({J/m3s})$',fontsize=fonts1)
# ax2.set_ylim(298, 330)

# plt.xlim(0.008, 0.01)
plt.tick_params(labelsize=10)

# handles1, label1 = ax1.get_legend_handles_labels()
# handles2, label2 = ax2.get_legend_handles_labels()
ax1.legend(loc='lower left', fontsize=12)
ax2.legend(loc='upper left', fontsize=12)

plt.show()