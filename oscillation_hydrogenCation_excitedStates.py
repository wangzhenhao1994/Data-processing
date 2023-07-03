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
omega = [588.26682,1346.61947,1479.61893,1588.32041,1704.69493,1821.06945,1945.11702,2053.8185,2189.37563,3294.29418,3530.87974,3771.30183,4161.34831,357.69977,815.04447,4518.73414,4718.93935,4911.44435,5110.36619,5463.29203,3081.36344]
amplitude=[0.5,1,3,7,9,8,3,1,2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,]
y=0
i=0
for w in omega:
    y = y+ np. np.cos(t/33356.40952*w*2*np.pi)*amplitude[i]
plt.plot(t,y)
plt.show()