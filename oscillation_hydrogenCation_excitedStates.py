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

t=np.arange(0,1300,0.001)
omega = [526.94738,811.57589,884.65618,1343.65171,1597.50957,1697.51418,2155.22759,2328.3125,2678.32864,3212.96868,3656.57888,1146.2067]
amplitude=[0.5,1,3,7,9,8,3,1,2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,]
y=0
i=0
for w in omega:
    y = y+ np. np.cos(t/33356.40952*w*2*np.pi)*amplitude[i]
plt.plot(t,y)
plt.show()