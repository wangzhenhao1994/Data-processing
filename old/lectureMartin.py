from turtle import color
import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import speed_of_light as lc
print(lc)

SMALL_SIZE = 15
MEDIUM_SIZE = 22
BIGGER_SIZE = 26

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

x=np.arange(-60,60,0.01)
sigma = 30/2/np.sqrt(2*np.log(2))
y=np.exp(-x*x/2/sigma/sigma)#*np.cos(x/2.6*np.pi*2)
dy=x/sigma/sigma*np.exp(-x*x/2/sigma/sigma)

fig = plt.figure()
gs = fig.add_gridspec(3, hspace=0)
axs = gs.subplots(sharex=True, sharey=True)
#fig.suptitle('Sharing both axes')
axs[0].plot(x, y, color='r')
axs[0].set(ylabel='Intensity')
axs[1].plot(x, -y+0.7, color='b')
axs[1].set(ylabel='Phase')
axs[2].plot(x, dy*10+0.3,color='b')
axs[2].set(ylabel='Frequency', xlabel=r'$\tau$ (fs)')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
plt.yticks([], [])
plt.show()