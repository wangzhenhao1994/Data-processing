from matplotlib.colors import Normalize
from matplotlib import cm

from calculate_k_b import Calibration_mass
import numpy as np
import scipy as sp
import scipy.signal as sps
import scipy.interpolate as spi
import obspy
from obspy.signal.detrend import polynomial
import matplotlib as mpl
from matplotlib import pyplot as plt
import h5py
import pathlib as pl
import os
import pickle

def save_obj(obj, filename):
    with open(filename, 'wb+') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

dataPath = r'C:\Users\user\Desktop\Data_newTOF'
yearAndMonth = r'2021-10'
date, filename = [r'2021-10-10', r'scan_tof_2021-10-10-22-28-56']
dataPath = pl.PureWindowsPath(dataPath)
if not os.path.exists(os.path.join(dataPath, yearAndMonth, date, filename)):
    os.mkdir(os.path.join(dataPath, yearAndMonth, date, filename))
dataFile =  os.path.join(dataPath, yearAndMonth, date, filename+r'.hdf5')
massDictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_massDict'+r'.pkl')
interSpecDictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_interSpec'+r'.pkl')
interSpecDictFile2 = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_interSpec2'+r'.pkl')
FFTdictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_FFTdict'+r'.pkl')
SFTdictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_SFTdict'+r'.pkl')
delayFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_delay'+r'.pkl')
delayDictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_delayDict'+r'.pkl')
files = [dataFile, massDictFile, interSpecDictFile, interSpecDictFile2, FFTdictFile, SFTdictFile, delayFile, delayDictFile]

delayMonitor = load_obj(delayDictFile)
print(delayMonitor[0].shape)
#plt.plot(np.arange(delelyMonitor[0].shape[0])/1000,delelyMonitor[0])
windowsize=7
vmax=abs(delayMonitor[0]).max()
vmin=abs(delayMonitor[0]).min()
levels = np.arange(vmin,vmax,(vmax-vmin)/500)
norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
fre, t, sft = sps.stft(delayMonitor[0], fs=1000, noverlap=windowsize-2, nperseg=windowsize, nfft=windowsize*3)
plt.contourf(t,fre*0.001*6.67/2.11,sft, cmap='jet')
plt.show()
