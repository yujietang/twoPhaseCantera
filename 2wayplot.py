
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

filename0 = './result/nosprayflame/cantera_nheptane_nospray_phi0.9.csv'
filename1 = './result/nosprayflame/Rochette_nosprayflame_phi0.9.csv'
filename2 = './result/lean_phio=0.9/2way_7.csv'
filename3 = './result/lean_phio=0.9/Rochette_d20_phi0.9.csv'
filename4 = './result/rich_phio=1.3/2way_46.csv'

filename5 = './result/rich_phio=1.3/Rochette_nosprayflame_phi1.3.csv'
filename6 = './result/rich_phio=1.3/Rochette_sprayflame_rich.csv'
filename7 = './result/rich_phio=1.3/cantera_nheptane_nospray_phi1.3.csv'

filename8 = './result/2way_53.csv'
filename9 = './result/2way_1.csv'

fonts1 = 18
fonts2 = 25
figs = (13,9)

fig = plt.figure(figsize=figs)
plt.rcParams['font.family'] = 'Serif'

from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'

data0 = np.loadtxt(filename0, delimiter=',', skiprows=1)
data0 = np.transpose(data0)
z0 = data0[0]
Tg0 = data0[1]

data1 = np.loadtxt(filename1, delimiter=',', skiprows=1)
data1 = np.transpose(data1)
z1 = data1[0]
Tg1 = data1[1]

data2 = np.loadtxt(filename2, delimiter=',', skiprows=1)
data2 = np.transpose(data2)
z2 = data2[0]
Tg2 = data2[3]
Sh2 = data2[-2]

data3 = np.loadtxt(filename3, delimiter=',', skiprows=1)
data3 = np.transpose(data3)
z3 = data3[0]
Tg3 = data3[1]

data4 = np.loadtxt(filename4, delimiter=',', skiprows=1)
data4 = np.transpose(data4)
z4 = data4[0]
Tg4 = data4[3]
Sh4 = data4[-2]

data5 = np.loadtxt(filename5, delimiter=',', skiprows=1)
data5 = np.transpose(data5)
z5 = data5[0]
Tg5 = data5[1]

data6 = np.loadtxt(filename6, delimiter=',', skiprows=1)
data6 = np.transpose(data6)
z6 = data6[0]
Tg6 = data6[1]

data7 = np.loadtxt(filename7, delimiter=',', skiprows=1)
data7 = np.transpose(data7)
z7 = data7[0]
Tg7 = data7[1]


ax1  = fig.add_subplot(111)
# ax1.scatter(xp, d2/d2[0], label='$d$',c ='r',s=5)
# ax1.scatter(xp, d, label='$d$',c ='r',s=5)
# ax1.scatter(z, Sm, label='${S_m}$', c='r', s=3)

####lean fuel
# ax1.plot(z2-0.0017, Tg2, label='$d=20$ ${\mu}m$', c='b',linewidth = 3)
# ax1.plot(z0-0.0017, Tg0, label='$no$ $spray$', c='r', linewidth=3)
# ax1.scatter(z3, Tg3, label='$Rochette:$ $d=20$ ${\mu}m$', c = 'b')
# ax1.scatter(z1, Tg1, label='$Rochette:$ $no$ $spray$', c = 'k')

###rich fuel
ax1.plot(z4-0.0017, Tg4, label='$d=20$ ${\mu}m$', c='b',linewidth = 3)
ax1.scatter(z5, Tg5, label='$Rochette:$ $no$ $spray$', c = 'k')
ax1.scatter(z6, Tg6, label='$Rochette:$ $d=20$ ${\mu}m$', c = 'b')
ax1.plot(z0-0.0017, Tg0, label='$no$ $spray$', c='r', linewidth=3)
ax1.set_xlabel(r'$z [m]$',fontsize=fonts1)
# ax1.set_ylabel(r'$d^2/{d_0}^2 (-)$',fontsize=fonts1)
ax1.set_ylabel(r'$T_{gas}$ $({K})$',fontsize=fonts1)

# ax1.set_ylim(-2e-5, 0.0)
# ax1.set_xlabel(r'$z$ $(m)$',fontsize=fonts1)
# ax1.set_ylabel(r'$d$ $({\mu}m)$',fontsize=fonts1)

ax2 = ax1.twinx()
# ax2.scatter(xp, Tp/Tp[0], label='Droplet Temperature',c ='b',s=5)
# ax2.scatter(xp, Tp, label='Droplet Temperature',c ='b',s=5)
# ax2.scatter(z, Sh, label='${S_h}$',c ='b',s=3)
# ax2.plot(z2-0.0017, Sh2, label='$energy$ $source$',c ='y', linewidth = 3)
ax2.plot(z4-0.0017, Sh4, label='$energy$ $source$',c ='y', linewidth = 3)
ax2.set_ylabel(r'${S_h}$ $({J/m^3s})$',fontsize=fonts1)
# ax2.set_ylim(298, 330)

plt.xlim(0.006, 0.015)
plt.tick_params(labelsize=10)

# handles1, label1 = ax1.get_legend_handles_labels()
# handles2, label2 = ax2.get_legend_handles_labels()
ax1.legend(loc='upper left', fontsize=12)
ax2.legend(loc='lower right', fontsize=12)

plt.show()