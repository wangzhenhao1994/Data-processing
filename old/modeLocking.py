import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import speed_of_light as lc
print(lc)

SMALL_SIZE = 20
MEDIUM_SIZE = 22
BIGGER_SIZE = 26

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

L=5000*1e-9
x=np.arange(0,L,L/300000)

omega=lambda n: 2*np.pi*lc/L*n
wave = lambda n,x,t: np.sin(2*np.pi*x/L*n)*np.cos(omega(n)*t)
ax = plt.subplot(111)

t = lambda n: 2*np.pi/omega(n)
stack = lambda n,t: sum([wave(i,x,t) for i in range(n-50,n+50)],0)

ax.plot(x,stack(10000,t(10000)*3000),linewidth=3.0)

#ax.plot(x,stack(10000,t(10000)*6.5),linewidth=3.0)

#ax.plot(x,stack(20,t(20)*12),linewidth=3.0)
#ax.plot(x,stack(20,t(20)*13),linewidth=3.0)
#ax.plot(x,stack(20,t(20)*14),linewidth=3.0)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
#plt.xlabel('t/fs')

plt.yticks([])
plt.xticks([])
plt.savefig('pulse',transparent = True)
plt.show()
