import numpy as np
import scipy as sp
import scipy.fft
from matplotlib import pyplot as plt

SMALL_SIZE = 18
MEDIUM_SIZE = 22
BIGGER_SIZE = 26

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def FFT(t, y):
    n = len(t)
    delta = (max(t) - min(t)) / (n-1)
    k = int(n/2)
    f = np.arange(k) / (n*delta) / 10**12  # frequency unit THz
    Y = abs(sp.fft.fft(y))[:k]
    return np.array([f, Y])

x=np.arange(0,100*np.pi,np.pi/30)

y=np.concatenate((np.sin(x[:1500]),np.sin(x[1500:])))

for d in np.arange(0,np.pi,np.pi/5):
    y2 = np.concatenate((np.zeros(5000),np.sin(x[:1500]),np.sin(x[1500:]+np.pi*d),np.zeros(5000)))
    #y2 = np.concatenate((np.sin(x[:1500]),np.sin(x[1500:]+np.pi*d)))
    t,fre=FFT(x,y2)
    plt.plot(t,fre,label=str(d/np.pi)+str('\u03C0'))
    #plt.plot(x,y2,label=str(d/np.pi)+str('\u03C0'))

plt.legend()
plt.title("Phase mismatch simulation")
plt.show()

