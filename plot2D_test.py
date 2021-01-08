
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

filename1 = './test_origin.csv'
filename2 = './test.csv'

fonts1 = 23
fonts2 = 20
linew = 3
figs = (13,9)

plt.figure(figsize=figs)
plt.rcParams['font.family'] = 'Serif'
# plt.rcParams['font.serif'] = 'Dejavu Serif'
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'

data1 = np.loadtxt(filename1, delimiter=',', skiprows=1)
data2 = np.loadtxt(filename2, delimiter=',', skiprows=1)

data1 = np.transpose(data1)
data2 = np.transpose(data2)

time = data1[0]
xp = data1[1]
d = data1[2]
d2 = data1[3]
T = data1[-2]
ug = data1[-1]
sc = plt.scatter(time, d2,label=r'$diameter = 25{\mu}m$, Euler', marker='x',c ='r',s=16)

time = data2[0]
xp = data2[1]
d = data2[2]
d2 = data2[3]
T = data2[-2]
ug = data2[-1]
sc = plt.scatter(time, d2,label=r'$diameter = 25{\mu}m$, RK4', marker='.',c ='b',s=1,alpha=1)

# plt.xlim(0,0.05)
# plt.ylim(0,2600)
plt.tick_params(labelsize=fonts2)
# plt.xlabel(r'$time [s]$',fontsize=fonts2)
plt.xlabel(r'$time$ $( s )$',fontsize=fonts2)
plt.ylabel(r'$d^2$ $( {\mu}m^2 )$',fontsize=fonts2)
plt.rcParams.update({'font.size':18})
plt.legend(loc='upper right')

plt.show()
