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

x=np.arange(-15,45,0.01)
sigma = 8/2/np.sqrt(2*np.log(2))
y=np.exp(-x*x/2/sigma/sigma)*np.cos(x/5*np.pi*2)+np.exp(-(x-30)*(x-30)/2/sigma/sigma)*np.cos((x-30)/5*np.pi*2)
ax = plt.subplot(111)
ax.plot(x,y,color = 'r',linewidth=10.0)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
#plt.xlabel('phase')

plt.yticks([])
plt.xticks([])
plt.tight_layout()
plt.savefig('pulse',transparent = True)
plt.show()