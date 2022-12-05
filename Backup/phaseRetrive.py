from operator import xor
from pickle import FALSE
#import matlab.engine
#import matlab
import matplotlib as mpl
from matplotlib import pyplot as plt
mpl.style.use('seaborn-poster')
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib import gridspec
import numpy as np
import scipy as sp
import scipy.signal as sps
import scipy.fft as sft
import scipy.interpolate as spi
import pywt
import matplotlib as mpl
from obspy.signal.detrend import polynomial
import h5py
import pathlib as pl
import os
from math import ceil, pi, sqrt, log 
import pickle
from brokenaxes import brokenaxes
from adjustText import adjust_text
from decimal import Decimal
from cal_intensity import cal_intensity
from calculate_k_b import Calibration_mass



t = np.linspace(8.96215610638882E-14, 1.2665985222923E-12, 2969)
omega = [524,809,879,1340,1593,1688,2149,2673,3207,3647]
amp=[1,1,1,1,1,1,1,1,1,10]
def phaseRetrive(omega,t,amp):
    osc=np.sum(np.dot(np.reshape(amp,(1,len(amp))),np.array([np.cos(2*np.pi*f*t) for f in np.array(omega)*29979245800])),0)
    return osc

plt.plot(t,phaseRetrive(omega,t,amp))
plt.show()