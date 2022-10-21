from fileinput import filename
from operator import xor
from pickle import FALSE
from certifi import where
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn import gaussian_process
from yaml import load
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
from Superlets.python.superlet import *
import time

import originpro as op
ifsave = False
ifsaveT = False
ifsavePhase = False
if ifsave or ifsaveT or ifsavePhase:
    op.set_show(show=True)
op.set_show(show=True)

my_path = os.path.abspath(__file__)


mpl.rcParams['lines.linewidth'] = 1
plt.rcParams["font.family"] = "arial"
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.it'] = 'Arial:italic'
mpl.rcParams['mathtext.rm'] = 'Arial'


SMALL_SIZE = 26
MEDIUM_SIZE = 22
BIGGER_SIZE = 26

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def save_obj(obj, filename):
    with open(filename, 'wb+') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

class FFT_ionS():
    def __init__(self,folder):
        '''

        '''
        self.filepath = []
        self.rootPath= pl.PureWindowsPath(r'C:\Users\user\Desktop\Data_newTOF\dataProcessing\CO2')
        #self.savePath= pl.PureWindowsPath(r'C:\Users\user\Desktop\Data_newTOF\dataProcessing\4.5E+14_H2O')
        self.savePath= pl.PureWindowsPath(os.path.join(r'C:\Users\user\Desktop\Data_newTOF\dataProcessing\CO2',folder))
        self.delayB = np.linspace(start=0,stop=1300, num=13000)*10**-15
        self.delayStep = 0.1*10**-15
        self.delayB_noPad = np.linspace(start=0,stop=1300, num=13000)*10**-15
        self.rebin_delay = None
        self.ifrebin = False
        self.longstage = True
        self.data = 0
        self.specBigBottleB = {}
        self.specBigBottleB_noPad = {}
        self.interSpectraBottleB = {}
        self.phaseBottleB = {}
        self.filterPhase = {}
        self.filterSpec = {}
        self.fftSB = {}
        self.stftSB = {}
        self.specBigBottleB = {}
        self.waveletsCoe = {}
        self.zeroIndex = {}
        self.zeroIndexMin = 10000000000000000000000000
        self.zeroIndexMax = 0
        self.retrivedPhase = {}
        self.spec = None
        self.scanNo = None
        self.stage = 'piezo'
        self.intensity = folder
        self.molecule = str('CO2')

        #setting for data from digitizer
        self.channelSize = 12032#24000#1536
        self.scanLengthD = 3320#1200
        self.peakRange = [-100, 100]  # range of the peak
        self.delayD = np.arange(self.scanLengthD)/self.scanLengthD*100*2*2*3.33564*10**-15*3602/3647
        self.spectraBottleD = {}
        self.fftSD = {}
        self.stftSD = {}
        self.dataD = 0

    def read(self):
        plt.clf()
        totalCounter = {}
        #for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
        #    self.fftSB[gas] = 0
        for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
            self.specBigBottleB[gas] = np.zeros((1000,13000))
            totalCounter[gas] = 0
        if not os.path.exists(os.path.join(self.savePath,r'plotData')):
            os.mkdir(os.path.join(self.savePath,r'plotData'))
        for file in os.listdir(self.savePath):
            filename = os.fsdecode(file)
            if filename.endswith(".pkl"):
                self.interSpectraBottleB = load_obj(os.path.abspath(os.path.join(self.savePath, filename)))
                for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
                    for i in range(self.interSpectraBottleB[gas].shape[0]):
                        self.specBigBottleB[gas][totalCounter[gas]] = self.interSpectraBottleB[gas][i]
                        totalCounter[gas] = totalCounter[gas]+1
        print(totalCounter)
        for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
            self.specBigBottleB[gas] = self.specBigBottleB[gas][:totalCounter[gas]]
            if self.specBigBottleB[gas].shape[0] == 1:
                self.specBigBottleB[gas] = self.specBigBottleB[gas][0].reshape(1,self.specBigBottleB[gas].shape[1])
                self.specBigBottleB[gas] = np.flip(self.specBigBottleB[gas],axis=1)
                #plt.plot(self.specBigBottleB[gas])
                #plt.show()
            else:
                zeroIndex2 = self.findZeroDelay2()
                zeroIndex2Min = np.amin(zeroIndex2)
                zeroIndex2Max = np.amax(zeroIndex2)
                for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
                    for i in range(self.specBigBottleB[gas].shape[0]):
                        self.specBigBottleB[gas][i] = np.concatenate((self.specBigBottleB[gas][i][(zeroIndex2[i]-zeroIndex2Min):], np.full(zeroIndex2[i]-zeroIndex2Min,self.specBigBottleB[gas][i][-1])))
                    
    def interFFT(self, t, y): #https://www.researchgate.net/post/What-are-the-basic-differences-between-FFT-and-DFT-and-DCT
        n = len(t)
        #y=np.roll(y,int(n/2)) #circular shift
        delta = (max(t) - min(t)) / (n)
        f = sft.fftshift(sft.fftfreq(n, delta))/ 10**12  # frequency unit THz
        fft_y = sft.fftshift(sft.fft(y))
        fft_y = fft_y#/np.max(np.abs(fft_y))
        #P = np.angle(fft_y)%(2*np.pi)
        Y = np.abs(fft_y)
        return np.array([f, Y])

    def FFT(self): #https://www.researchgate.net/post/What-are-the-basic-differences-between-FFT-and-DFT-and-DCT# do FFT for total spectra
        for gas in self.specBigBottleB.keys():
            #_,y,t = self.inter_window(self.specBigBottleB[gas],self.delayB_noPad,windowSize=200)
            y,t = self.inter_padding(self.specBigBottleB[gas],self.delayB,paddingSize=100000)
            self.fftSB[gas] = self.interFFT(t, y)           
    #def FFT2(self): #https://www.researchgate.net/post/What-are-the-basic-differences-between-FFT-and-DFT-and-DCT # do FFT for seperate time trace and add together.
    #    for gas in self.interSpectraBottleB.keys():
    #        for i in range(self.interSpectraBottleB[gas].shape[0]):
    #            _,y,t = self.inter_window(self.interSpectraBottleB[gas][i],self.delayB_noPad,windowSize=200)
    #            y,t = self.inter_padding(y,t,paddingSize=100000)
    #            f,Y = self.interFFT(t, y)
    #            self.fftSB[gas] = self.fftSB[gas] + Y

    def transition(self):
        self.specBigBottleB_noPad=self.specBigBottleB.copy()
        #self.specBigBottleB_noPad['Ch2'] = self.specBigBottleB_noPad['Ch2']-self.specBigBottleB_noPad['Ch10']
        self.delayB_noPad = self.delayB.copy()
        if not os.path.exists(os.path.join(self.savePath,r'totalSpec.pkl')):
            save_obj(self.specBigBottleB,os.path.join(self.savePath,r'totalSpec.pkl'))
            save_obj(self.delayB_noPad,os.path.join(self.savePath,r'delay.pkl'))
        else:
            print('Total spectra Saved!\n')
        for gas in self.specBigBottleB.keys():
            self.specBigBottleB[gas] = self.specBigBottleB[gas].sum(axis = 0)#np.take(self.specBigBottleB[gas],[i for i in range(self.specBigBottleB[gas].shape[0]) if i not in self.goodSpecIndex],axis=0).sum(axis = 0)
        self.specBigBottleB['Ch2'] = self.specBigBottleB['Ch2']-self.specBigBottleB['Ch0']
        self.specBigBottleB['Ch4'] = self.specBigBottleB['Ch4']-self.specBigBottleB['Ch10']
        #self.specBigBottleB['Ch6'] = self.specBigBottleB['Ch6']-self.specBigBottleB['Ch10']
        #self.specBigBottleB['Ch8'] = self.specBigBottleB['Ch8']-self.specBigBottleB['Ch10']

    def checkSavedData(self):
        if os.path.exists(os.path.join(self.savePath,r'totalSpec.pkl')):
            self.specBigBottleB = load_obj(os.path.join(self.savePath,r'totalSpec.pkl'))
            self.delayB = load_obj(os.path.join(self.savePath,r'delay.pkl'))
            print('Saved total spectra detected!\n')
            try:
                self.phaseBottle = load_obj(os.path.join(self.savePath, r'phaseBottle.pkl'))
                self.interPhase = load_obj(os.path.join(self.savePath, r'interPhase.pkl'))
                print('Saved phase detected!\n')
            except:
                print('No saved phase detected!\n')
            return 0
        else:
            return 1
    
    def window(self, windowSize=0, direction='left'):
        '''
        windowSize is in fs.
        '''
        windowSize = windowSize+np.abs(self.delayB[0])*1e15
        self.windowSize = windowSize
        for gas in list(self.specBigBottleB):
            if 'rebin' not in gas:
                self.specBigBottleB[gas], self.specBigBottleB['window_'+gas], self.delayB = self.inter_window(self.specBigBottleB[gas], self.delayB, windowSize=windowSize, direction=direction)
            else:
                self.specBigBottleB[gas], self.specBigBottleB['window_rebin_'+gas], self.rebin_delay = self.inter_window(self.specBigBottleB[gas], self.rebin_delay, windowSize=windowSize, direction=direction)
        for gas in list(self.spectraBottleD):
            if 'rebin' not in gas:
                self.spectraBottleD[gas], self.spectraBottleD['window_'+gas], self.delayD = self.inter_window(self.spectraBottleD[gas], self.delayD, windowSize=windowSize, direction=direction)
            else:
                self.spectraBottleD[gas], self.spectraBottleD['window_rebin_'+gas], self.rebin_delay = self.inter_window(self.spectraBottleD[gas], self.rebin_delay, windowSize=windowSize, direction=direction)
            
    def inter_window(self, data, delay, windowSize=0, direction='left', useWindow = True):
        '''
        windowSize is in fs.
        '''
        __len = np.size(data)
        windowSize = int(windowSize*1e-15/(delay[1]-delay[0]))
        if direction == 'left':
            window = data[-(__len-windowSize):]
            delay = delay[-(__len-windowSize):]
        elif direction == 'right':
            window = data[:__len-windowSize]
            delay = delay[:__len-windowSize]
        elif direction == 'middle':
            window = np.concatenate((
                data[:round((__len-windowSize)/2)],  # need to check again
                (data[round((__len-windowSize)/2)] + \
                 data[round((__len+windowSize)/2)])/2*(np.zeros(windowSize)+1),
                data[-(__len-round((__len+windowSize)/2)+1):]
            ))
        window=polynomial(window, order=15, plot=False)
        hwindow = np.hanning(len(window))
        if useWindow:
            window2=window*hwindow #https://stackoverflow.com/questions/55654699/how-to-get-correct-phase-values-using-np-fft
        else:
            window2=window
        return window, window2, delay

    def padding(self, paddingSize=0):
        for gas in list(self.specBigBottleB):
            if "rebin" not in gas:
                self.specBigBottleB[gas], self.delayB = self.inter_padding(self.specBigBottleB[gas], self.delayB, paddingSize = paddingSize)
            else:
                self.specBigBottleB[gas], self.rebin_delay = self.inter_padding(self.specBigBottleB[gas], self.rebin_delay, paddingSize = paddingSize)
        for gas in list(self.spectraBottleD):
            if "rebin" not in gas:
                self.spectraBottleD[gas], self.delayD = self.inter_padding(self.spectraBottleD[gas], self.delayD, paddingSize = paddingSize)
            else:
                self.spectraBottleD[gas], self.rebin_delay = self.inter_padding(self.specBigBottleB[gas], self.rebin_delay, paddingSize = paddingSize)

    def inter_padding(self, inter_data, delay, paddingSize=0):
        data = None
        delayStep = delay[90]-delay[89]
        delay = np.concatenate((
            #np.arange(delay[0]-(paddingSize+1) *
            #          delayStep, delay[0], delayStep),
            delay,
            np.arange(delay[-1]+delayStep, delay[-1] +
                      (paddingSize)*delayStep, delayStep)
        ))
        
        data = np.concatenate(
            #(np.zeros(paddingSize)+inter_data[0], inter_data, np.zeros(paddingSize)+inter_data[-1]), axis=0)
            (inter_data, (np.zeros(paddingSize))), axis=0)
        delay = delay[:len(data)]
        return data, delay

    def interpS(self, interNum=1):
        interpDelay = np.arange(
            self.delayB[0], self.delayB[-1], (self.delayB[1]-self.delayB[0])/(interNum+1))
        for gas in self.specBigBottleB.keys():
            iS = np.interp(interpDelay, self.delayB, self.specBigBottleB[gas])
            self.specBigBottleB[gas] = iS
            # print(np.shape(iS))
        self.delayB = interpDelay

    def FFTS(self):
        for gas in self.specBigBottleB.keys():
            if "rebin" not in gas:
                delay = self.delayB
                self.fftSB[gas] = self.FFT(
                    delay, self.specBigBottleB[gas])
            else:
                delay = self.rebin_delay
                self.fftSB[gas] = self.FFT(
                    delay, self.specBigBottleB[gas])
        for gas in self.spectraBottleD.keys():
            if self.ifrebin:
                if "rebin" not in gas:
                    continue
                else:
                    delay = self.rebin_delay
                    self.fftSD[gas] = self.FFT(
                        delay, self.spectraBottleD[gas])
            else:
                if "rebin" not in gas:
                    delay = self.delayD
                    self.fftSD[gas] = self.FFT(
                        delay, self.spectraBottleD[gas])
                else:
                    continue

    def rmvExp(self):
        for gas in self.specBigBottleB.keys():
            y = self.specBigBottleB[gas]
            a = np.array([[self.delayB[0], 1], [self.delayB[-1], 1]])
            b = np.array([y[0],y[-1]])
            [k, b] = np.linalg.solve(a, b)
            self.specBigBottleB[gas] = self.specBigBottleB[gas]-(k*self.delayB+b)
            #self.specBigBottleB[gas]=polynomial(self.specBigBottleB[gas], order=3, plot=False)
            hwindow = np.hanning(len(y))#sps.windows.flattop(len(window),sym=FALSE)
            self.specBigBottleB[gas]=self.specBigBottleB[gas]*hwindow

    def runingAverage(self, n=5):
        def runAve(x, N): return np.convolve(x, np.ones(N)/N, mode='valid')
        for gas in self.specBigBottleB.keys():
            self.specBigBottleB[gas] = runAve(self.specBigBottleB[gas], n)
            new_size = len(self.specBigBottleB[gas])
        self.delayB = self.delayB[:new_size]

    def smooth(self, windowSize=100, order=3):
        for gas in self.specBigBottleB.keys():    
            self.specBigBottleB[gas]=sps.savgol_filter(self.specBigBottleB[gas], windowSize, order) # window size 51, polynomial order 3
    
    def smooth2(self,data,windowSize=100,order=3):
        return sps.savgol_filter(data, windowSize, order)

    def rebin_factor(self, a, factor):
            '''Rebin an array to a new shape.
            newshape must be a factor of a.shape.
            '''
            newshape = tuple(int(i/factor) for i in a.shape)
            i1 = tuple(int(i*factor) for i in newshape)
            a = a[:i1[0]]
            return a.reshape(newshape+(-1,)).mean(axis=-1)
    
    def rebinS(self, factor=10):
        self.ifrebin = True
        for gas in list(self.specBigBottleB.keys()):
            self.specBigBottleB['rebin_'+gas] = self.rebin_factor(self.specBigBottleB[gas], factor)
        self.rebin_delay = self.rebin_factor(self.delayB, factor)

    def dataProcess(self):
        #self.mass_spectra()
        self.fftSB()

    def useFilter(self, lowcut, highcut):
        for gas in list(self.specBigBottleB.keys()):
            # if 'rebin' not in gas:
            fs = 1/(self.delayB[1]-self.delayB[0])
            # else:
            #     fs = 1/(self.rebin_delay[90]-self.rebin_delay[89])
            self.specBigBottleB['filter_'+gas] = self.butter_bandpass_filter(self.specBigBottleB[gas], lowcut, highcut, fs)

    def useFilter2(self, data, lowcut, highcut):
        # if 'rebin' not in gas:
        fs = 1/(self.delayB[1]-self.delayB[0])
        # else:
        #     fs = 1/(self.rebin_delay[90]-self.rebin_delay[89])
        return self.butter_bandpass_filter(data, lowcut, highcut, fs)

    def butter_bandpass_filter(self, data, lowcut, highcut, fs, order=2):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = sps.butter(order, [low, high], analog=False, btype='band', output='sos')
        y = sps.sosfiltfilt(sos, data)
        #w, h = sps.sosfreqz(sos, worN=2000)
        #plt.plot(w, h)
        #plt.show()
        return y

    def findZeroDelay(self):
        zeroIndex = []
        for i in range(self.interSpectraBottleB['Ch8'].shape[0]):
            _Spec = self.interSpectraBottleB['Ch8'][i]
            _Spec = self.butter_bandpass_filter(_Spec, 500/33.35641*1e12, 2000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            #plt.plot(_Spec)
            #plt.show()
            zeroIndex =zeroIndex + [(_Spec[:1500]).argmax()]
        return zeroIndex

    def findZeroDelay2(self):
        zeroIndex = []
        for i in range(self.specBigBottleB['Ch8'].shape[0]):
            _Spec = self.specBigBottleB['Ch8'][i]
            #_Spec = self.specBigBottleB['Ch8']
            _Spec = self.butter_bandpass_filter(_Spec, 500/33.35641*1e12, 2000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            plt.plot(_Spec[:2000])
            plt.show()
            zeroIndex =zeroIndex + [(_Spec[:2000]).argmax()]
        return zeroIndex

    def findZeroDelay3(self,j):
        _Spec = 0
        for gas in self.interSpectraBottleB.keys():
            _Spec = _Spec+self.interSpectraBottleB[gas][j]
        _Spec = self.butter_bandpass_filter(_Spec, (12493-3000)/33.35641*1e12, (12493+3000)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        #plt.plot(np.abs(_Spec[:1500]))
        #plt.show()
        zeroIndex =(np.abs(_Spec[:1500])).argmax()
        return zeroIndex

    def findZeroDelay4(self):
        _Spec = self.specBigBottleB['Ch8']
        _Spec = self.butter_bandpass_filter(_Spec, 10/33.35641*1e12, 2000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        #plt.plot(_Spec[:1000])
        #plt.plot(_Spec[-1000:])
        #plt.show()
        #zeroIndex =np.abs(_Spec[-1000:]).argmax()
        zeroIndex =np.abs(_Spec[:1000]).argmax()
        return zeroIndex

    def show_FFT(self):
        self.dcRange=0
        gs = gridspec.GridSpec(2, 3)
        fig = plt.figure(figsize=(15,5))
        lab=['Mass 16','Mass 22','Mass 28','Mass 44','ref16','ref22']
        i=0
        if ifsave:
            self.wks = op.new_sheet('w',lname=str('FFT')+str('_')+self.molecule+str('_')+self.intensity)
        for gas in ['Ch2','Ch4','Ch6','Ch8','Ch0','Ch10']:
            #if 'filter' in gas  or 'window' in gas or 'rebin' in gas or 'com' in gas:
            #   continue
            #elif gas not in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch0+6']:
            #    continue

            ax = fig.add_subplot(gs[i])
            axP = ax.twinx()
            label=lab[i]
            #[f, Y, _] = self.fftSB[gas]
            #[f_filter, Y_filter, P_filter] = self.fftSB['filter_'+gas]
            [f_window, Y_window] = self.fftSB['window_'+gas]
            #[f_rebin, Y_rebin, P_rebin] = self.fftSB['rebin_'+gas]
            plotRatio=1
            f_window=f_window.astype(float)*33.35641
            #f_rebin=f_rebin*33.35641
            #Y_window=np.where(f_window>=100,Y_window,0)/np.amax(Y_window)
            
            f_window = f_window[(f_window>100)]
            
            aa=f_window.shape[0]
            bb=Y_window.shape[0]
            Y_window = Y_window[bb-aa:bb]
            Y_window=Y_window/np.amax(Y_window)


            #P_window = np.sum([np.where(np.logical_and(f_window>f-20,f_window<f+20),P_window,0)
            # for f in [525.62772,805.21693,883.50191,1338.30036,1591.79458,1696.17455,2147.24515,3205.9563,3645.84333,2669.14501]],0)
            #P_window = np.where(Y_window>0.2,P_window,0)
            #Y_filter=np.where(np.logical_and(f>=10, f<=4500),Y_filter,0)/plotRatio
            #ax.plot(f,np.abs(Y), label=label)#gas)#, label=self.trueScanTime)
            #ax.plot(f,np.abs(Y_filter)/20, label='filter'+label)#gas)#, label=self.trueScanTime)
            ax.plot(f_window, np.abs(Y_window), label=label)#gas)#, label=self.trueScanTime)
            #ax.plot(f_rebin,np.abs(Y_rebin), label='rebin_window_'+label)#gas)#, label=self.trueScanTime)
            #ax.plot(f,np.real(Y), label="Re")#gas)#, label=self.trueScanTime)
            #ax.plot(f,np.imag(Y), label="Im")#gas)#, label=self.trueScanTime)
            #axP.plot(f_window,P_window/np.pi,'r')

            #ax.plot(f[peaksIndex[i]],np.abs(Y_window)[peaksIndex[i]], "rx")
            #if gas == 'Ch2+4':
            #    [f_window, Y_ref, P_window] = self.fftSB['window_'+'Ch2-4']
            #    ax.plot(f, np.abs(Y_ref), label='Ch2-4')
            ax.set_xlabel('Frequency/cm-1')
            ax.set_ylabel('Amplitude')
            #axP.set_ylabel('Phase/pi')
            ax.legend(loc=2)
            #ax.set_ylim([0,15])

            ax.set_xlim([450,14000])
            #ax.set_xlim([58870,59920])
            #ax.set_title(gas)
            i=i+1
            if ifsave:
                self.wks.from_list(0, f_window, lname="Frequency", axis='X')
                self.wks.from_list(i, np.abs(Y_window), lname=label, axis='Y')

        plt.legend()
        #fig.tight_layout()
        plt.show()

    def show_FFTD(self):
        #plt.figure()
        #unused_f, ref = self.fftS['CH3OH+']
        i=12
        for gas in self.fftSD.keys():
            #if 'filter' not in gas:
            #    continue
            [f, Y, _] = self.fftSD[gas]
            #Y=np.where(np.logical_and(f>=10, f<=4500),Y,0)
            if gas == 'window_H2+':
                plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=r'$\mathrm{H_2^+}$')
            #if gas == 'H+':
            #    plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=r'$\mathrm{H^+}$')
            if gas == 'window_O+':
                plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=r'$\mathrm{O^+}$')
            if gas == 'window_O+H2':
                plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=r'$\mathrm{H_2^++O^+}$')
                print("The max of the FFT is ",np.amax( Y[self.dcRange:]/1000))
            if gas == 'window_O-H2':
                plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=r'$\mathrm{H_2^+-O^+}$')
            else:
                #plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=gas)
                pass
            plt.xlabel(r'$\mathrm{Frequency (cm^{-1})}$')
            plt.ylabel('a.u.')
            plt.legend(loc=2)
            #plt.xlabel('Frequency/cm-1')
            #plt.ylim([0,np.max(Y[self.dcRange:1000])*3/2])
            plt.xlim([0,4500])
            i=i+1
            if ifsave:
                self.wks.from_list(7, f[self.dcRange:], 'X')
                self.wks.from_list(i, np.abs(Y[self.dcRange:]), lname=gas, axis='Y')
        #plt.legend()
        plt.tight_layout()
        d='_100_120'

    def show_Spectra(self, shift=0):
        gs = gridspec.GridSpec(2, 3)
        #gs = gridspec.GridSpec(1, 1)
        fig = plt.figure(figsize=(20,8))
        #ax = fig.add_subplot(111)
        lab=['Mass 16','Mass 22','Mass 28','Mass 44']
        i=0
        if ifsaveT:
            self.wks = op.new_sheet('w',lname=str('Time')+str('_')+self.molecule+str('_')+self.intensity)
        delay = (self.delayB-self.findZeroDelay4()*self.delayStep)
        for gas in ['Ch2','Ch4','Ch6','Ch8']:
            
            #if 'filter' in gas or 'window' in gas or 'rebin' in gas or 'com' in gas:
            #    continue
            #if gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
            #    pass
            #else:
            #    print('Skip!')
            #    continue
            label=lab[i]
            
            ax = fig.add_subplot(gs[i])
            i=i+1
            
            rebindelay = self.rebin_delay
            #else:
            #    delay = self.rebin_delay
            ax.plot(delay*10**15+shift,
                     self.specBigBottleB[gas]-np.mean(self.specBigBottleB[gas][:9000]), label=gas)
            #ax.plot(delay*10**15+shift,
            #         self.specBigBottleB['filter_'+gas], label=gas)#/np.amax(self.specBigBottleB[gas])
            #ax.plot(delay*10**15+shift,
            #         self.specBigBottleB['window_filter_'+gas], label=label)
            #ax.plot(delay*10**15+shift,
            #         self.specBigBottleB['window_'+gas], label='window_'+gas)
            #ax.plot(rebindelay*10**15+shift,
            #         self.specBigBottleB['rebin_window_'+gas], label='rebin_window_'+gas)
            ax.set_xlabel("Delay/fs")
            ax.set_ylabel('a.u.')
            #plt.xlim([300,400])
            ax.legend()
            print(gas)
            if ifsaveT:
                self.wks.from_list(0, delay*10**15, lname='Delay (fs)',axis='X')
                self.wks.from_list(i, (self.specBigBottleB[gas]-np.mean(self.specBigBottleB[gas][-9000:]))/(np.amax(self.specBigBottleB[gas][-9000:])-np.amin(self.specBigBottleB[gas][-9000:])), lname=label, axis='Y')
        #plt.legend()
        fig.tight_layout()
        plt.show()

    def calDrift(self, _cRange, gas='Ch8'):
        plt.clf()
        print(self.specBigBottleB[gas].shape)
        _iS = np.array(np.zeros(self.specBigBottleB[gas].shape))
        for i in range(self.specBigBottleB[gas].shape[0]):
            #_iS[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], (3643.52883-0.1)/33.35641*1e12, (3643.52883+0.1)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            _iS[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], 500/33.35641*1e12, 2000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        try:
            iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
            _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=20000)
            indexMax = []
            for i in range(self.specBigBottleB[gas].shape[0]):
                _inter=iS.ev(i,_delayRange)
                plt.plot(_delayRange,_inter)
                #indexMax = indexMax + [sps.argrelextrema(_inter, np.greater)[0][-1]]
                indexMax = indexMax + [np.argmax(np.abs(_inter))]
            plt.xlabel('Delay (fs)')
            plt.ylabel('a.u.')
            plt.show()

            _ref = indexMax[int(self.specBigBottleB[gas].shape[0]/2)]
            _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/20000
            for i in range(self.interSpectraBottleB[gas].shape[0]):
                _inter=iS.ev(i,_delayRange+_shift[i])
                plt.plot(_delayRange,_inter)
            plt.xlabel('Delay (fs)')
            plt.ylabel('a.u.')
            plt.show()
        except:
            plt.clf()
            empty = np.zeros((3,_iS.shape[1]))
            _iS = np.concatenate((_iS, empty), axis=0)
            iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
            _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=20000)
            indexMax = []
            for i in range(self.specBigBottleB[gas].shape[0]):
                _inter=iS.ev(i,_delayRange)
                plt.plot(_delayRange,_inter)
                #indexMax = indexMax + [sps.argrelextrema(_inter, np.greater)[0][-1]]
                indexMax = indexMax + [np.argmax(_inter)]
            plt.xlabel('Delay (fs)')
            plt.ylabel('a.u.')
            plt.show()
    
            _ref = indexMax[int(self.specBigBottleB[gas].shape[0]/2)]
            _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/20000
            for i in range(self.interSpectraBottleB[gas].shape[0]):
                _inter=iS.ev(i,_delayRange+_shift[i])
                plt.plot(_delayRange,_inter)
            plt.xlabel('Delay (fs)')
            plt.ylabel('a.u.')
            plt.show()

        return _shift

    def delayCorrection(self, _cRange = [0,150]):
        xxx = 0
        while True:
            try:
                _shift = self.calDrift(_cRange = _cRange)
                break
            except IndexError:
                print('Wrong Region! Try again!')
                if xxx>10:
                    _cRange = [np.array(_cRange)[0], np.array(_cRange)[1]+0.5]
                _cRange = np.array(_cRange)+0.5
                xxx = xxx+1
        
        iinter = {}
        oShape = self.specBigBottleB['Ch6'].shape[0]
        for gas in self.specBigBottleB.keys():
            if _shift.size ==1:
                empty = np.zeros((3,self.specBigBottleB[gas].shape[1]))
                self.specBigBottleB[gas] = np.concatenate((self.specBigBottleB[gas], empty), axis=0)
            else:
                continue
            iS = sp.interpolate.RectBivariateSpline(range(self.specBigBottleB[gas].shape[0]), self.delayB_noPad*1e15, self.specBigBottleB[gas])
            inter = np.zeros((oShape,self.specBigBottleB[gas].shape[1]))
            for i in range(oShape):
                inter[i] = iS.ev(i, self.delayB_noPad*1e15+_shift[i])
            iinter[gas] = inter
        self.specBigBottleB = iinter

    def phaseRetrive(self,omega=[3643.52883]):
        #[520.50412,619.22042,804.6874,879.47248,1136.73313,1334.16573,1588.43498,1687.15128,2144.83594,2318.33731,2668.33146,3200.80119,3643.52883]
        self.interSpectraBottleB = {}
        for gas in self.specBigBottleB_noPad.keys():
            _shape = self.specBigBottleB_noPad[gas].shape[0]
            _shape = _shape-_shape%6
            self.interSpectraBottleB[gas] = self.specBigBottleB_noPad[gas][:_shape]
            self.interSpectraBottleB[gas] = self.interSpectraBottleB[gas].reshape(int(_shape/6), 6, self.interSpectraBottleB[gas].shape[1]).sum(axis=1)
            
            self.filterPhase[gas]={}
            self.filterSpec[gas]={}
            for i in range(self.interSpectraBottleB[gas].shape[0]):
                self.filterPhase[gas][str(i)]={}
                self.filterSpec[gas][str(i)]={}
                for f in omega:
                    f0=(f-1)/33.35641*1e12
                    f1=(f+1)/33.35641*1e12
                    #dd,t=self.inter_padding(self.interSpectraBottleB[gas][i],self.delayB,paddingSize=10000000)
                    self.filterSpec[gas][str(i)][str(f)] = self.useFilter2(self.interSpectraBottleB[gas][i],f0,f1)
                    #if gas=='Ch6' and f==3643.52883:
                    #    plt.plot(self.filterSpec[gas][str(i)][str(f)]*3000-382,label='noPad')
                    #    dd,t=self.inter_padding(self.interSpectraBottleB[gas][i],self.delayB,paddingSize=10000000)
                    #    plt.plot(self.useFilter2(dd,f0,f1)-382,label='Pad')
                    #    plt.plot(self.interSpectraBottleB[gas][i])
                    #    plt.legend()
                    #    plt.show()
                    self.filterPhase[gas][str(i)][str(f)] = np.angle(sps.hilbert(self.filterSpec[gas][str(i)][str(f)]))
        #zeroIndex = self.findZeroDelay()
        #sampleIndex = np.zeros(self.delayB_noPad.size)+1
        #sampleIndex = np.where(np.abs(self.filterPhase['Ch8'][str(i)][str(f)]-np.pi)<2*np.pi*2/9.15*(self.delayB[1]-self.delayB[0])*1e15,sampleIndex,0)
        #sampleIndex = np.where(self.delayB_noPad*1e15<900, sampleIndex, 0)
        #sampleIndex = np.where(self.delayB_noPad*1e15>200, sampleIndex, 0)
        for gas in self.interSpectraBottleB.keys():
            self.phaseBottleB[gas] = {}
            for f in omega:
                self.phaseBottleB[gas][str(f)] = {}
                for i in range(self.interSpectraBottleB[gas].shape[0]):
                    zeroIndex = self.findZeroDelay3(i)
                    print(zeroIndex)
                    #delay = (self.delayB_noPad-zeroIndex[i]*self.delayStep)*10**15
                    delay = (self.delayB_noPad-zeroIndex*self.delayStep)*10**15
                    delay = delay[(delay>200)]
                    phi = np.unwrap(self.filterPhase[gas][str(i)][str(f)])
                    phi = phi[-delay.size:]
                    delay = delay[(delay<200+9.15497*60)]
                    phi = phi[:delay.size]
                    #plt.plot(phi)
                    #plt.plot(delay/(33356.40952/f)*2*np.pi)
                    #plt.show()
                    phi = phi-delay/(33356.40952/f)*2*np.pi
                    #plt.plot(phi)
                    #plt.show()
                    phi = phi%(2*np.pi)
                    if i==0:
                        inter = np.zeros((self.interSpectraBottleB[gas].shape[0],phi.size))
                        inter[0] = phi
                    else:
                        inter[i] = phi
                self.phaseBottleB[gas][str(f)]=inter
                    #if phi[]
        for gas in self.phaseBottleB.keys():
            if gas == 'Ch10':
                continue
            for f in omega:
                plt.plot(np.mean(self.phaseBottleB[gas][str(f)],axis=1), label = str(f))
            plt.legend()
            plt.show()
        #save_obj(self.)
    
    def fit_sin(self, tt, yy, fre):
        '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
        tt = np.array(tt)
        yy = np.array(yy)
        ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
        Fyy = abs(np.fft.fft(yy))
        guess_freq = abs(ff[np.argmax(Fyy[2:])+2])   # excluding the zero frequency "peak", which is related to offset
        guess_amp = np.std(yy) * 2.**0.5
        guess_offset = np.mean(yy)
        guess = np.array([guess_amp, 2.*np.pi*fre, 0., guess_offset])

        def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
        popt, pcov = sp.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
        A, w, p, c = popt
        f = w/(2.*np.pi)
        return popt

    def phaseRetrive2(self,omega=[1283.3, 1379]):
        """
        Retrieve the absolute phase
        """
        self.interSpectraBottleB = {}
        self.phaseBottle = {}

        for gas in ['Ch8','Ch2','Ch4','Ch6']:
            _shape = self.specBigBottleB_noPad[gas].shape[0]
            if _shape == 1:
                self.interSpectraBottleB[gas] = self.specBigBottleB_noPad[gas][:_shape]
                _zeroIndex = [self.findZeroDelay4()]
            else:
                _shape = _shape-_shape%6
                self.interSpectraBottleB[gas] = self.specBigBottleB_noPad[gas][:_shape]
                self.interSpectraBottleB[gas] = self.interSpectraBottleB[gas].reshape(int(_shape/6), 6, self.interSpectraBottleB[gas].shape[1]).sum(axis=1)
                _zeroIndex = self.findZeroDelay()
            self.phaseBottle[gas] = {}
            for j in range(self.interSpectraBottleB[gas].shape[0]):
                self.phaseBottle[gas][str(j)]={}
                interPhase=[]
                
                _,_ss,_d = self.inter_window(self.interSpectraBottleB['Ch8'][j],self.delayB_noPad,150)
                y,t = self.inter_padding( _ss,_d,paddingSize=100000)
                [f_window, Y_window] = self.interFFT(t, y)
                f_window = f_window*33.35641
                Y_window = Y_window[:len(f_window)]
                amp = np.array([np.mean(np.where(np.logical_and(f_window>=ff-20, f_window<=ff+20),Y_window,0)) for ff in omega])
                amp = np.reshape(amp,(1,len(amp)))

                _,_s,_d = self.inter_window(self.interSpectraBottleB[gas][j],self.delayB_noPad,150)
                y,t = self.inter_padding( _s,_d,paddingSize=100000)
                [f_window, Y_window] = self.interFFT(t, y)
                f_window = f_window*33.35641
                osc = _s
                hwindow = np.hanning(len(_d))
                cal_len = int(33356.40952/1000/self.delayStep/1e15*2)#choose how much steps to calculate
                for i in range(cal_len):
                    step=self.delayB[1]-self.delayB[0]
                    oscA=np.sum(np.dot(amp,np.array([np.sin(2*np.pi*f*(_d-(_zeroIndex[j]+i)*step)) for f in np.array(omega)/33.35641*1e12])),0)#
                    #sps.windows.flattop(len(window),sym=FALSE)
                    #plt.plot(self.delayB_noPad,hwindow*(oscA/100))
                    #plt.plot(self.delayB_noPad,osc)
                    #plt.show()
                    y,t=self.inter_padding(hwindow*(oscA/np.amax(oscA))+osc/np.amax(polynomial(osc, order=15, plot=False)), _d, paddingSize = 100000)
                    f_inter,Y = self.interFFT(t,y)
                    #plt.plot(f_inter*33.35641,Y)
                    #plt.xlim([300,4500])
                    #plt.show()
                    interPhase = interPhase+[Y]
                
                #plt.contourf(interPhase)
                #plt.show()
                f_inter = f_inter*33.35641
                index0 = np.argmin(np.abs(f_inter-300))
                index1 = np.argmin(np.abs(f_inter-4500))
                self.wks = op.new_sheet('w',lname=str('Axis'))
                self.wks.from_list(0, f_inter[np.argmin(np.abs(f_inter-10)):np.argmin(np.abs(f_inter-15000))], lname='Frequency (cm-1)', axis='Y')
                self.wks.from_list(1, range(cal_len)*self.delayStep*1e15, lname='Time (fs)', axis='Y')
                self.interPhase = np.array(interPhase)[:,np.argmin(np.abs(f_inter-10)):np.argmin(np.abs(f_inter-15000))]
                save_obj(self.interPhase, os.path.join(self.savePath, r'interPhase.pkl'))
                #if ifsavePhase:
                #    self.wks = op.new_sheet('w',lname=str('Phase')+str('_')+str(gas)+str('_')+str('%.1E' % Decimal(self.intensity))+str('_'))
                for f in omega:
                    self.phaseBottle[gas][str(j)][str(f)] = np.sum(np.array(interPhase)[:,np.argmin(np.abs(f_inter-f+5)):np.argmin(np.abs(f_inter-f-5))],1)
                    tt = self.delayB[:cal_len]-_zeroIndex[j]*self.delayStep
                    res = self.fit_sin(tt,self.phaseBottle[gas][str(j)][str(f)],f/33.35641*1e12)
                    self.phaseBottle[gas][str(j)][str(f)+'_fit'] = [tt, res]
                    A, w, p, c = self.phaseBottle[gas][str(j)][str(f)+'_fit'][1]
                    fitfunc = lambda t: A * np.sin(w*t + p) + c
                    if f == 3643.52883:
                        plt.plot(tt, fitfunc(tt), "r-", label="y fit curve", linewidth=2)
                        plt.plot(tt,self.phaseBottle[gas][str(j)][str(f)])
            #plt.show()
        save_obj(self.phaseBottle, os.path.join(self.savePath, r'phaseBottle.pkl'))

    def superlet(self,gas):
        for gas in gas:
            #fs = 1/(self.delayB[1]-self.delayB[0])
            fs = 1/(self.rebin_delay[1]-self.rebin_delay[0])  # sampling frequency
            signal = self.specBigBottleB['rebin_'+'window_'+gas]
            # frequencies of interest in Hz
            foi = np.linspace(int(300/33.35641)*1e12, int(4500/33.35641)*1e12, 300)
            scales = scale_from_period(1 / foi)
            spec = superlet(
                signal,
                samplerate=fs,
                scales=scales,
                order_max=100,
                order_min=1,
                c_1=0.3,
                adaptive=True,
            )
            # amplitude scalogram
            ampls = np.abs(spec)

            fig, (ax1, ax2) = plt.subplots(2, 1,
                                           sharex=True,
                                           gridspec_kw={"height_ratios": [1, 3]},
                                           figsize=(6, 6))

            ax1.plot(np.arange(signal.size) / fs, signal, c='cornflowerblue')
            ax1.set_ylabel('signal (a.u.)')

            extent = [0, len(signal) / fs, foi[0], foi[-1]]
            im = ax2.imshow(ampls, cmap="jet", aspect="auto", extent=extent, origin='lower')

            plt.colorbar(im,ax = ax2, orientation='horizontal',
                         shrink=0.7, pad=0.2, label='amplitude (a.u.)')

            ax2.plot([0, len(signal) / fs], [20, 20], "--", c='0.5')
            ax2.plot([0, len(signal) / fs], [40, 40], "--", c='0.5')
            ax2.plot([0, len(signal) / fs], [60, 60], "--", c='0.5')

            ax2.set_xlabel("time (s)")    
            ax2.set_ylabel("frequency (Hz)")

            fig.tight_layout()
            plt.show()

    def wlt(self, gas):
        PHztocm = 33356.40952
        filt_freqs_in_cm = np.array([3643.52883])
        bw = 30
        wavelet = 'cmor3-1'
        signal = self.specBigBottleB['rebin_'+'window_'+gas] - np.mean(self.specBigBottleB['rebin_'+'window_'+gas])
        signal = signal/np.max(signal)
        delay = self.rebin_delay*1e15
        scales = np.logspace(1, 2.7, num=400)

        frequencies = pywt.scale2frequency(
            wavelet, scales) / (delay[1] - delay[0])
        print('min freq (cm-1): ' + str(frequencies.min()*PHztocm))
        print('max freq (cm-1): ' + str(frequencies.max()*PHztocm))

        [coefficients, frequencies] = pywt.cwt(
            signal, scales, wavelet, delay[1] - delay[0])
        amplitude = np.abs(coefficients)
        amplitude = amplitude/np.max(amplitude)

        fig, ax = plt.subplots(
            2+filt_freqs_in_cm.size,
            1,
            gridspec_kw={'height_ratios': [3, 1, ] + filt_freqs_in_cm.size*[1]})

        im = ax[0].pcolormesh(
            delay,
            frequencies*PHztocm,
            #10*np.log10(amplitude),
            amplitude,
            cmap='jet',
            vmax = 0.5)
        ax[0].set_xlim([delay.min(), delay.max()])
        ax[0].set_ylim([0, 7000])
        ax[0].set_ylabel('frequency [cm-1]')
        ax[0].xaxis.set_ticklabels([])
        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes('top', size='7%', pad='2%')
        cb = fig.colorbar(im, cax=cax, orientation='horizontal')
        cb.set_label('amplitude [dB]')
        cax.xaxis.set_ticks_position('top')
        cax.xaxis.set_label_position('top')
        ax[0].plot([0], [0], label=str(gas))
        ax[0].legend(frameon=0, loc='upper right')

        norm = 0
        for idx, filt_freq_in_cm in enumerate(filt_freqs_in_cm):
            filter = (
                (frequencies*PHztocm > filt_freq_in_cm - bw/2)
                & (frequencies*PHztocm < filt_freq_in_cm + bw/2))
            occ = amplitude[filter, :].sum(axis=0)
            norm = np.max([occ.max(), norm])

        for idx, filt_freq_in_cm in enumerate(filt_freqs_in_cm):
            filter = (
                (frequencies*PHztocm > filt_freq_in_cm - bw/2)
                & (frequencies*PHztocm < filt_freq_in_cm + bw/2))
            occ = amplitude[filter, :].sum(axis=0)
            occ = occ/norm
            ax[1+idx].plot(
                delay,
                10*np.log10(occ),
                label='ampl. [dB], freq = ' + str(filt_freq_in_cm))
            ax[1+idx].set_xlim([delay.min(), delay.max()])
            ax[1+idx].set_ylim([-30, 1])
            ax[1+idx].xaxis.set_ticklabels([])
            ax[1+idx].legend(frameon=0, loc='upper right')

        ax[-1].plot(delay, signal)
        ax[-1].set_xlim([delay.min(), delay.max()])
        ax[-1].set_ylim([
            np.mean(signal)-2*np.std(signal),
            np.mean(signal)+2*np.std(signal)])
        ax[-1].set_ylabel('mass amplitude')
        ax[-1].set_xlabel('time (fs)')

        # fig.tight_layout()
        #plt.savefig('wavelet_transform_mass_' + str(gas) + '.png')
        plt.show()

    def plotPhase(self):
        if ifsavePhase:
            self.wks = op.new_sheet('w',lname=str('Phase')+str('_')+str('H2O')+str('%.1E' % Decimal(self.intensity)))
        self.phase = {}
        ii=1
        for gas in ['Ch2','Ch4','Ch6','Ch8']:#
            self.phase[gas] = {}
            #if gas == 'Ch8':
            #    omega = [520.50412,619.22042,804.6874,879.47248,1136.73313,1334.16573,1588.43498,1687.15128,2144.83594,2318.33731,2668.33146,3200.80119,3643.52883]
            #elif gas == 'Ch0':
            #    omega = [520.50412,804.6874,1334.16573,1588.43498,2144.83594,3643.52883]
            #elif gas == 'Ch2' or gas == 'Ch4':
            #    omega = [3643.52883]
            #elif gas == 'Ch6':
            #    omega = [520.50412,804.6874,1334.16573,1588.43498,2144.83594,3643.52883]
            omega = [1283.3, 1379]
            mean = []
            std = []
            for f in omega:#520.50412,619.22042,804.6874,879.47248,1136.73313,1334.16573,1588.43498,1687.15128,2144.83594,2318.33731,2668.33146,3200.80119,
                _inter = []
                for i in self.phaseBottle[gas].keys():
                    #if f==3643.52883:
                    #    plt.plot(self.phaseBottle[gas][str(i)][str(f)])
                    #    plt.show()
                    #A, w, p, c = self.phaseBottle[gas][str(i)][str(f)+'_fit'][1]
                    p = sps.argrelextrema(self.phaseBottle[gas][str(i)][str(f)], np.greater)[0][0]/(33356.40952/f/self.delayStep/1e15)*np.pi*2
                    _inter = _inter + [p]
                #plt.plot(_inter)
                #plt.show()
                _std = np.std(_inter)
                _mean = np.mean(_inter)
                mean = mean + [_mean]
                std = std + [_std]
            self.phase[gas] = [omega, mean, std]
            if ifsavePhase:
                self.wks.from_list(0, omega, 'X')
                self.wks.from_list(ii, mean, lname=gas, axis='Y')
                self.wks.from_list(ii+1, std, lname=gas, axis='E')
                ii=ii+2
            plt.errorbar(omega,np.array(mean)/np.pi,yerr=np.array(std)/np.pi, fmt='s',label = gas)
        plt.legend()
        plt.show()

    def phaseRetrive3(self):
        """
        Retrieve the relative phase to the phase in 'Ch8'
        """
        omega=[520.50412,619.22042,804.6874,879.47248,1136.73313,1334.16573,1588.43498,1687.15128,2144.83594,2318.33731,2668.33146,3200.80119,3643.52883]
        self.interSpectraBottleB = {}
        self.phaseBottle = {}

        for gas in ['Ch8','Ch0','Ch2','Ch4','Ch6']:
            _shape = self.specBigBottleB_noPad[gas].shape[0]
            _shape = _shape-_shape%6
            self.interSpectraBottleB[gas] = self.specBigBottleB_noPad[gas][:_shape]
            self.interSpectraBottleB[gas] = self.interSpectraBottleB[gas].reshape(int(_shape/6), 6, self.interSpectraBottleB[gas].shape[1]).sum(axis=1)
            print(self.findZeroDelay())
            self.phaseBottle[gas] = {}
            for j in range(1,self.interSpectraBottleB[gas].shape[0]):
                print(gas)
                self.phaseBottle[gas][str(j)]={}
                interPhase=[]
                
                _,_ss,_d = self.inter_window(self.interSpectraBottleB['Ch8'][j], self.delayB_noPad,150, useWindow = False)
                _,_s,_d = self.inter_window(self.interSpectraBottleB[gas][j], self.delayB_noPad,150, useWindow = False)

                cal_len = int(33356.40952/520.50412/self.delayStep/1e15*2)#choose how much steps to calculate
                osc = _s[:(len(_s)-cal_len)]/np.amax(_s)
                hwindow = np.hanning(len(osc))
                
                for i in range(cal_len):
                    oscA=_ss[i:(len(_ss)-cal_len+i)]/np.amax(_ss)*0.2
                    y,t=self.inter_padding(hwindow*(osc+oscA), _d, paddingSize = 100000)
                    f_inter,Y = self.interFFT(t,y)
                    #plt.plot(f_inter*33.35641,Y)
                    #plt.xlim([300,4500])
                    #plt.show()
                    interPhase = interPhase+[Y]
                f_inter = f_inter*33.35641
                for f in omega:
                    self.phaseBottle[gas][str(j)][str(f)] = np.sum(np.array(interPhase)[:,np.argmin(np.abs(f_inter-f+5)):np.argmin(np.abs(f_inter-f-5))],1)
                    tt = np.array(range(cal_len))*self.delayStep
                    res = self.fit_sin(tt,self.phaseBottle[gas][str(j)][str(f)],f/33.35641*1e12)
                    self.phaseBottle[gas][str(j)][str(f)+'_fit'] = [tt, res]
                    A, w, p, c = self.phaseBottle[gas][str(j)][str(f)+'_fit'][1]
                    fitfunc = lambda t: A * np.sin(w*t + p) + c
                    #plt.plot(tt, fitfunc(tt), "r-", label="y fit curve", linewidth=2)
                    #plt.plot(tt,self.phaseBottle[gas][str(j)][str(f)])
                    #plt.show()
        save_obj(self.phaseBottle, os.path.join(self.savePath, r'phaseBottle.pkl'))

    def plotInterPhase(self):
        self.wks = op.new_sheet('w',lname=str('FFT')+str('_')+self.molecule+str('_')+self.intensity)
        #print(len(self.interPhase),len(self.interPhase[0]))
        self.wks.from_list2(self.interPhase)

if __name__ == '__main__':
    #for ff in [r'pu1.2E+15pr7.4E+13_CO2', r'pu2.5E+14pr7.4E+13_CO2', r'pu4.7E+14pr7.4E+13_CO2', r'pu8.4E+14pr7.4E+13_CO2']:
    for ff in [r'pu7.4E+13pr4.7E+14_CO2', r'pu7.4E+13pr8.4E+14_CO2', r'pu7.4E+13pr1.2E+15_CO2']:
    #for ff in [r'pu4.7E+14pr1.2E+15_CO2', r'pu4.7E+14pr8.4E+14_CO2', r'pu4.7E+14pr6.5E+14_CO2', r'pu4.7E+14pr4.2E+14_CO2', r'pu4.7E+14pr2.6E+14_CO2']:
        d = FFT_ionS(ff)#pu1.2E+15pr7.4E+13_CO2,pu2.5E+14pr7.4E+13_CO2,pu4.7E+14pr7.4E+13_CO2,pu8.4E+14pr7.4E+13_CO2
        if d.checkSavedData():
            d.read()
            d.delayCorrection()
        d.transition()
        #d.delayB = np.flip(d.delayB)
        #for gas in d.specBigBottleB.keys():
        #    d.specBigBottleB[gas] = np.flip(d.specBigBottleB[gas])
            #plt.plot(d.delayB,d.specBigBottleB[gas])
            #plt.show()
        
        d.phaseRetrive2()
        #d.plotPhase()
        #for i in range(d.specBigBottleB_noPad['Ch8'].shape[0]):
        #    plt.plot(d.specBigBottleB_noPad['Ch8'][i])
        #plt.show()
        #d.smooth(9)
        d.plotInterPhase()
        #d.window(windowSize=150, direction='left')
        #d.useFilter(10/33.35641*1e12, 6000/33.35641*1e12)
        #d.show_Spectra()
        #d.FFT()

        #d.show_FFT()
        #d.rebinS(5)
        #d.superlet(['Ch8'])
        #d.wlt('Ch6')
