import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
fonts1 = 18
fonts2 = 25
figs = (13,9)

fig = plt.figure(figsize=figs)
plt.rcParams['font.family'] = 'Serif'

plt.rcParams['mathtext.fontset'] = 'stix'

for n in np.arange(1, 8, 1):
    filename0 = './result/1way_'+str(n)+'.csv'
    # filename1 = './LPT1.csv'
    data1 = np.loadtxt(filename0, delimiter=',', skiprows=1)
    data1 = np.transpose(data1)
    #2-way:
    # z = data1[0]
    # ug = data1[1]
    # rhog = data1[2]
    # Tg = data1[3]
    # Yf = data1[4]
    # YO2 = data1[5]
    # YOH = data1[6]
    # YCO = data1[7]
    # YCO2 = data1[8]
    # phig = data1[-4]
    # Sm = data1[-3]
    # Sh = data1[-2]
    # Syf = data1[-1]

    #1-way:
    time = data1[0]
    xp = data1[1]
    d = data1[2]
    d2 = data1[3]
    mp = data1[4]
    Tp = data1[5]
    T = data1[6]
    mtfd = data1[-6]
    htfd = data1[-5]
    Sm = data1[-4]
    Sh = data1[-3]
    rhop = data1[-2]
    Mliquid = data1[-1]

    # plt.plot(z, Sm, label='${Y_f}$'+str(n))
    plt.scatter(n, Mliquid[0], label='${liquid mass flux}$'+str(n))
    # plt.plot(xp, Sm, label='${S_m}-$'+str(n))


plt.xlabel(r'$z [m]$',fontsize=fonts1)
plt.ylabel(r'${Y_f}$ $({-})$',fontsize=fonts1)
plt.tick_params(labelsize=10)

plt.legend(loc='lower left', ncol=4, fontsize=12)

plt.show() 