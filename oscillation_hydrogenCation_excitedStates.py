import numpy as np
from matplotlib import pyplot as plt

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

t=np.arange(0,15,0.001)
omega = 3656.57888
#amplitude=[1]
i=0

y = 1.5+np.sin(t/33356.40952*omega*2*np.pi)#*amplitude[i]
y2 = 1.5+np.sin(t/33356.40952*omega*2*np.pi+np.pi/2)
y3 = 1.5+np.sin(t/33356.40952*omega*2*np.pi-np.pi/2)
plt.plot(t,y,'r',label=r'$\phi=0$')
plt.plot(t,y2,'b',label=r'$\phi=\pi/2$')
plt.plot(t,y3,'g',label=r'$\phi=-\pi/2$')
plt.title(r'Ion yield = Sin($\omega t+\phi$)')
plt.legend()
plt.show()