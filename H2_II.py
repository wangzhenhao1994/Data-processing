from cProfile import label
from fileinput import filename
from weakref import ref
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
import pycwt
import matplotlib as mpl
from obspy.signal.detrend import polynomial
import h5py
import pathlib as pl
import os
import pickle
from decimal import Decimal
from cal_intensity import cal_intensity
from calculate_k_b import Calibration_mass
from BaselineRemoval import BaselineRemoval
from lmfit import minimize, Parameters, Parameter, report_fit
import originpro as op

import math

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
        self.folder = folder
        self.filename = filename
        self.saveRef = str(filename)
        self.filepath = []
        self.delayB,self.stepSize = np.linspace(start=0,stop=1300, num=13000,endpoint=False,retstep=True)
        self.delayB = self.delayB*10**-15/4161.16632*4155.0587
        self.stepSize = self.stepSize*10**-15
        self.delayB_noPad = np.linspace(start=0,stop=1300, num=13000)*10**-15/4161.16632*4155.0587
        self.rebin_delay = None
        self.ifrebin = False
        self.longstage = True
        self.data = 0
        self.specBigBottleB = {}
        self.specBigBottleB_noPad = {}
        self.interSpectraBottleB = {}
        self.phaseSpecBottleB = {}
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
        self.intensity = 8.9E14
        self.molecule = str('H2')
        self.rootPath= pl.PureWindowsPath(r'D:\SF_FFT')
        self.savePath= pl.PureWindowsPath(os.path.join(r'D:\SF_FFT\processedData\inter',self.molecule,folder))
    def read(self):
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
        zeroIndex2 = self.findZeroDelay2()
        zeroIndex2Min = np.amin(zeroIndex2)
        zeroIndex2Max = np.amax(zeroIndex2)
        for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
            for i in range(self.specBigBottleB[gas].shape[0]):
                self.specBigBottleB[gas][i] = np.concatenate((self.specBigBottleB[gas][i][(zeroIndex2[i]-zeroIndex2Min):], np.full(zeroIndex2[i]-zeroIndex2Min,self.specBigBottleB[gas][i][-1])))
                    
    def interFFT(self, y):
        n = len(y)
        delta = self.stepSize#/(self.interNum+1)
        self.dw = 1/((n)*self.stepSize)/1E12*33.35641
        f = np.fft.fftfreq(n, delta)/1E12*33.35641  # frequency unit cm-1
        f = f[f>0]
        y=np.fft.fftshift(y)
        fft_y = np.fft.fft(y)
        #fft_y = fft_y*np.exp(1j*2*np.pi/n*int(3500+self.windowSize)*np.arange(n))#5415,
        fft_y=fft_y[:np.size(f)]
        #plt.plot(f,np.abs(fft_y))
        #plt.show()
        return f, fft_y

    def FFT3(self, windowSize = 100, delayRange=False, rebinF=1, paddingF=0, useWindow = False, zeroDirection = 'left', phaseCompensate = True, smooth=True, test=False):
        '''
        windowSize in fs, rebinF = inputShape/outputShape, padddingF means the length of the padding divided by length of the data.
        '''
        self.smoothedT={}
        self.smooth = smooth
        self.stepSize = self.stepSize*rebinF
        _size = self.specBigBottleB['Ch0'].size
        paddingSize = int(_size*paddingF)
        self.paddingSize = paddingSize
        [1,0.01,0.01,1,5,0.01]
        for gas in self.phaseSpecBottleB.keys():
            self.smoothedT[gas]=0
            for i in range(self.phaseSpecBottleB[gas].shape[0]):
                interSpec = self.phaseSpecBottleB[gas][i][-self.specBigBottleB[gas].size:]
                _,interSmoothedT,self.smoothedT['T'] = self.inter_window(interSpec,self.delayB,windowSize=50,direction=zeroDirection, useWindow=False, phaseCompensate=False,smooth = False)
                self.smoothedT[gas]= self.smoothedT[gas] + interSmoothedT
                if delayRange:
                    y=interSpec[np.where(np.logical_and(delayRange[0]<self.delayB,self.delayB<delayRange[1]))]
                    t=self.delayB[np.where(np.logical_and(delayRange[0]<self.delayB,self.delayB<delayRange[1]))]
                    _,y,t = self.inter_window(y,t,windowSize=0,direction=zeroDirection, useWindow=useWindow, phaseCompensate=phaseCompensate,smooth = smooth)
                else:
                    _,y,t = self.inter_window(interSpec,self.delayB,windowSize=windowSize,direction=zeroDirection, useWindow=useWindow, phaseCompensate=phaseCompensate,smooth = smooth)


                #_,y,t = self.inter_window(self.phaseSpecBottleB[gas][i],self.delayB,windowSize=windowSize,direction=zeroDirection, useWindow=useWindow, phaseCompensate=phaseCompensate,smooth = smooth)
                #plt.plot(t,y)
                #plt.plot(self.delayB,self.specBigBottleB['Ch8']-np.mean(self.specBigBottleB['Ch8']))
                #plt.show()
                y,t = self.inter_padding(y,t,paddingSize=paddingSize)
                if test:
                    y=np.sin(2*np.pi/300*np.arange(np.size(t)))+np.sin(2*np.pi/100*np.arange(np.size(t))+np.pi)+np.sin(2*np.pi/200*np.arange(np.size(t)))
                    y,t = self.inter_padding(y,t,paddingSize=paddingSize)
                #plt.plot(t,y)
                #plt.show()
                if rebinF < 1.5:
                    pass
                else:
                    y = self.rebin_factor(y,rebinF)
                    t = self.rebin_factor(t,rebinF)
                #self.interNum = 100
                #y = self.interInterp(t,y,self.interNum)
                f, Y = self.interFFT(y)
                if i == 0:
                    _interY = np.zeros((self.phaseSpecBottleB[gas].shape[0],Y.size),dtype=np.complex64)
                _interY[i] = Y
            self.fftSB['window_'+gas+'_fft'] = _interY
            self.fftSB[gas] = 0

        self.fftSB['window_'+'fre'] = f
        print(self.folder)
        save_obj(self.fftSB,os.path.join(self.savePath,str(smooth)+str(self.folder)+r'fftSB.pkl'))
        
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
            n=10
            self.phaseSpecBottleB[gas] = self.specBigBottleB[gas][-int(self.specBigBottleB[gas].shape[0]/n)*n:].reshape(int(self.specBigBottleB[gas].shape[0]/n), n, self.specBigBottleB[gas].shape[1]).sum(axis=1)
            self.specBigBottleB[gas] = self.specBigBottleB[gas].sum(axis = 0)#np.take(self.specBigBottleB[gas],[i for i in range(self.specBigBottleB[gas].shape[0]) if i not in self.goodSpecIndex],axis=0).sum(axis = 0)
        self.specBigBottleB['Ch2'] = self.specBigBottleB['Ch2']-self.specBigBottleB['Ch10']
        self.label = {}
        self.label['Ch0'] = 'Mass1'
        self.label['Ch2'] = 'Mass2'
        self.label['Ch4'] = 'Mass16'
        self.label['Ch6'] = 'Mass17'
        self.label['Ch8'] = 'Mass18'

    def checkSavedData(self):
        if os.path.exists(os.path.join(self.savePath,r'totalSpec.pkl')):
            self.specBigBottleB = load_obj(os.path.join(self.savePath,r'totalSpec.pkl'))
            self.delayB = load_obj(os.path.join(self.savePath,r'delay.pkl'))
            print('Saved total spectra detected!\n')
            try:
                self.phaseBottle = load_obj(os.path.join(self.savePath, r'phaseBottle.pkl'))
                print('Saved phase detected!\n')
            except:
                print('No saved phase detected!\n')
            return 0
        else:
            return 1
    
    def _flip(self):
        for gas in list(self.specBigBottleB):
            self.specBigBottleB[gas] = np.flip(self.specBigBottleB[gas],axis=0)

    def window(self, windowSize=0, direction='left', useWindow = True):
        '''
        windowSize is in fs.
        '''
        #windowSize = windowSize+np.abs(self.delayB[0])*1e15 #????????????????????????
        for gas in list(self.specBigBottleB):
            if 'rebin' not in gas:
                self.specBigBottleB[gas], self.specBigBottleB['window_'+gas], self.delayB = self.inter_window(self.specBigBottleB[gas], self.delayB, windowSize=windowSize, direction=direction, useWindow=useWindow)
            else:
                self.specBigBottleB[gas], self.specBigBottleB['window_'+gas], self.rebin_delay = self.inter_window(self.specBigBottleB[gas], self.rebin_delay, windowSize=windowSize, direction=direction, useWindow=useWindow)
    
    def window2(self, windowSize=0, direction='left'):
        '''
        windowSize is in fs.
        '''
        for gas in list(self.specBigBottleB):
            if 'rebin' not in gas:
                self.specBigBottleB[gas], self.specBigBottleB['window_'+gas], self.delayB = self.inter_window2(self.specBigBottleB[gas], self.delayB, windowSize=windowSize, direction=direction)
            else:
                self.specBigBottleB[gas], self.specBigBottleB['window_rebin_'+gas], self.rebin_delay = self.inter_window2(self.specBigBottleB[gas], self.rebin_delay, windowSize=windowSize, direction=direction)
        
    def inter_window(self, data, delay, windowSize=0, direction='left', useWindow = True, phaseCompensate = True, smooth = True):
        '''
        windowSize is in fs.
        '''
        data = data[:int(np.size(data)/2)*2]
        __len = np.size(data)
        windowSize = int(windowSize*1e-15/self.stepSize)
        self.windowSize = windowSize
        if smooth:
            data=sps.savgol_filter(data, window_length=10, polyorder=1,deriv=1,delta=10, mode='nearest')
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
        #window=polynomial(window, order=1, plot=False)
        window = window - window[0]#shift the base line to zero

        #if not math.log2(np.size(window)).is_integer():
        #    print('Fill the data to length if power of 2!')
        #    window = np.append(window,np.zeros(2**(math.ceil(math.log2(np.size(window))))-np.size(window)),axis=0)
        if useWindow:
            #window2=self.apply_triwindow(window)
            window2=self.apply_hannwindow(window)
            #window2=self.apply_hammwindow(window)
        else:
            window2=window


        if phaseCompensate:
            window2 = np.append(np.zeros(windowSize),window2,axis=0)
            window = np.append(np.zeros(windowSize),window,axis=0)

        return window, window2, delay

    def apply_hannwindow(self,y):
        hwindow = np.hanning(len(y))
        return y*hwindow

    def apply_hammwindow(self,y):
        hwindow = np.hamming(len(y))
        return y*hwindow

    def apply_triwindow(self,y):
        twindow = np.bartlett(len(y))
        return y*twindow

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
        delayStep = self.stepSize
        delay = np.concatenate((
            np.arange(delay[0]-(paddingSize+1) *
                      delayStep, delay[0], delayStep),
            delay,
            np.arange(delay[-1]+delayStep, delay[-1] +
                      (paddingSize)*delayStep, delayStep)
        ))
        data = np.concatenate(
            (np.zeros(paddingSize),inter_data, np.zeros(paddingSize)), axis=0)
            #(inter_data, np.zeros(paddingSize)), axis=0)
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

    def interInterp(self, t, S, interNum=1):
        interpDelay = np.arange(
            t[0], t[-1], self.stepSize/(interNum+1))
        return np.interp(interpDelay, self.delayB, S)
    def baseLineRemove(self,y):
        baseObj=BaselineRemoval(y)
        return baseObj.ZhangFit(lambda_=50, repitition=50)

    def rmvExp(self):
        for gas in self.specBigBottleB.keys():
            self.specBigBottleB[gas] = self.interRmvLinear(self.specBigBottleB[gas])

    def interRmvLinear(self,y):
        a = np.array([[0, 1], [np.size(y), 1]])
        b = np.array([y[0],y[-1]])
        [k, b] = np.linalg.solve(a, b)
        return y-(k*np.arange(np.size(y))+b)

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
            '''Rebin an array to a new shape.^^
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
        self.stepSize = self.stepSize*factor
        

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
        for i in range(self.interSpectraBottleB['Ch0'].shape[0]):
            _Spec = self.interSpectraBottleB['Ch0'][i]
            _Spec = self.butter_bandpass_filter(_Spec, 100/33.35641*1e12, 4000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            #plt.plot(_Spec)
            zeroIndex =zeroIndex + [(_Spec[:2000]).argmax()]
        #plt.show()
        return zeroIndex

    def findZeroDelay2(self):
        zeroIndex = []
        for i in range(self.specBigBottleB['Ch0'].shape[0]):
            _Spec = self.specBigBottleB['Ch0'][i]
            _Spec = self.butter_bandpass_filter(_Spec, 100/33.35641*1e12, 4000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            #plt.plot(_Spec)
            
            zeroIndex =zeroIndex + [np.abs((_Spec[:2000])).argmax()]
        #plt.show()
        return zeroIndex

    def findZeroDelay3(self):
        zeroIndex = []
        _Spec = self.specBigBottleB['Ch8']
        _Spec = self.butter_bandpass_filter(_Spec, 100/33.35641*1e12, 4000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        zeroIndex =np.abs((_Spec[:2000])).argmax()
        self.zeroIndex = zeroIndex # this number will be used for phase correction
        __inter = self.specBigBottleB_noPad.copy()
        self.specBigBottleB_noPad = {}
        for gas in self.specBigBottleB.keys():
            self.specBigBottleB[gas] = self.specBigBottleB[gas][zeroIndex:]
            self.specBigBottleB_noPad[gas]={}
            for j in range(__inter[gas].shape[0]):
                self.specBigBottleB_noPad[gas][j]=__inter[gas][j][zeroIndex:]
            #self.specBigBottleB[gas] = self.specBigBottleB[gas][:zeroIndex]
        self.delayB = np.arange(self.specBigBottleB['Ch8'].size)*self.stepSize
        return zeroIndex
    
    def phaseCorrection(self,spectra):
        spectra = self.interRmvLinear(spectra)
        interSpectra =spectra[:(self.zeroIndex*2)]
        L=np.size(spectra)-(self.zeroIndex*2)
        padded=np.concatenate((interSpectra[self.zeroIndex:self.zeroIndex*2],np.zeros(L),interSpectra[:self.zeroIndex]))
        Y=sp.fft.fft(padded)
        #plt.plot(np.unwrap(np.arctan(np.imag(Y)/np.real(Y))))
        #plt.show()
        F = sp.fft.idct(np.cos(np.unwrap(np.arctan(np.imag(Y)/np.real(Y)))),type=2,norm=None)+sp.fft.idst(np.sin(np.unwrap(np.arctan(np.imag(Y)/np.real(Y)))),type=2,norm=None)
        b=np.convolve(F,spectra)
        plt.plot(b)
        plt.show()
        #P_inter = np.unwrap(np.where(P_inter<0,P_inter+np.pi,P_inter))
        return b

    def show_FFT(self, ifsaveFFT=False, ifsaveT=False):
        self.dcRange=0
        fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True,figsize=(10,10))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel("Frequency ($\mathrm{cm^{-1}}$)")
        plt.ylabel("Amplitude (a.u.)", labelpad=25)
        plt.text(1.18,0.35,'Relative Phase v.s. Mass 2 ($\pi$)', rotation = 270)
        axPP = []
        for axx in ax.flatten():
            axPP = axPP + [axx.twinx()]
        axPP[0].get_shared_y_axes().join(axPP[0], axPP[1],axPP[2])

        lab=['Mass 2','Mass 1','Mass 18']
        i=0
        if ifsaveFFT:
            op.set_show(show=True)
            self.wksFFT = op.new_sheet('w',lname=str('FFT')+str('_')+self.folder)
        if ifsaveT:
            op.set_show(show=True)
            self.wksT = op.new_sheet('w',lname=str('Time')+str('_')+self.folder)

        _preP = 'window'+'_'
        self.result = {}
        self.result['phase'] = {}
        self.result['amplitude'] = {}

        boxcarA={}
        boxcarA['Ch0'],boxcarA['Ch2'],boxcarA['Ch4'],boxcarA['Ch6'],boxcarA['Ch8'] = [1,0.01,0.01,1,5]
        for gas in ['Ch2','Ch0','Ch8']:
            if 'filter' in gas  or 'window' in gas or 'rebin' in gas or 'com' in gas:
                continue
            elif gas not in ['Ch2','Ch0','Ch8']:
                continue
            
            axF=ax[i]
            axP = axPP[i]
            #axP.set_ylim([-np.pi,np.pi])
            label=lab[i]
            f_window = self.fftSB[_preP+'fre']
            Y = np.array([np.mean(self.fftSB[_preP+gas+'_fft'],axis=0), np.std(self.fftSB[_preP+gas+'_fft'],axis=0)])
            Y_window = np.array([np.mean(np.abs(self.fftSB[_preP+gas+'_fft']),axis=0),     np.std(np.abs(self.fftSB[_preP+gas+'_fft']),axis=0)])
            Y_window=np.where(f_window>=100,Y_window,0)/np.amax(Y_window)
            Y_window_im = np.array([np.mean(np.imag(self.fftSB[_preP+gas+'_fft']),axis=0),   np.std(np.imag(self.fftSB[_preP+gas+'_fft']),axis=0)])
            Y_window_re = np.array([np.mean(np.real(self.fftSB[_preP+gas+'_fft']),axis=0),   np.std(np.real(self.fftSB[_preP+gas+'_fft']),axis=0)])

            P_inter = np.angle(self.fftSB[_preP+gas+'_fft'])
            if gas == 'Ch2':
                P_ref = P_inter
                P_window =  np.array([np.mean(P_inter,axis=0), np.std(P_inter,axis=0)])
            else:
                P_inter = P_inter-P_ref
                for n in range(np.shape(P_inter)[0]):
                    P_inter[n]=np.where(P_inter[n]>2*np.pi-2, P_inter[n]-2*np.pi, P_inter[n])

                #P_window =  np.array([np.mean(P_inter,axis=0), np.sqrt(np.std(P_inter,axis=0)**2+np.std(P_ref,axis=0)**2)])
                P_window =  np.array([np.mean(P_inter,axis=0), np.sqrt(np.std(P_inter,axis=0)**2+0)])

            
            #aa = len(f_window[(f_window<100)])
            #bb=len(f_window[(f_window<4000)])
            aa = len(f_window[(f_window<100)])
            bb=len(f_window[(f_window<5000)])
            f_window = f_window[aa:bb]
            Y = Y[:,aa:bb]
            Y_window = Y_window[:,aa:bb]
            Y_window_im = Y_window_im[:,aa:bb]
            Y_window_re = Y_window_re[:,aa:bb]
            P_window = P_window[:,aa:bb]
            P_window[0]=np.where(P_window[0]<-np.pi-1,P_window[0]+2*np.pi,P_window[0])
            P_window[0]=np.where(P_window[0]>np.pi-1,P_window[0]-2*np.pi,P_window[0])
            P_window = P_window/np.pi
            self.result['frequency'] = f_window
            self.result['phase'][gas] = P_window
            omega = [1376.54213,1475.23277,1586.74037,1703.37476,1818.72746,1937.92524,2053.27793,2846.64812,3297.80532,3418.2848,3524.66561,3759.21609,3830.9911,3900.20271,4160.38712,4430.8251]
            for om in omega:
                if i>0:
                    axF.axvline(x=om,ymin=0,ymax=1.3,clip_on=False,c='k',linestyle='--',alpha=0.3)
                else:
                    axF.axvline(x=om,clip_on=False,c='k',linestyle='--',alpha=0.3)
                    axF.text(x=om-70, y=1.1, s=str(int(om)),fontsize='8')
                    if om==1594.43209:
                        axF.text(x=om-70, y=1.8, s='H2O+ Bend',fontsize='8')
                    if om==3655.52723:
                        axF.text(x=om-70, y=1.8, s='H2O a1 Sym',fontsize='8')
            
            inter = 0
            for w in omega:
                inter = inter + np.where(np.abs(f_window-w)<80,P_window,0)
            inter = np.where(inter==0,np.inf,inter)
            P_window = inter
            #fitRes = self.fit_phase(f_window,Y[0]/np.amax(np.abs(Y[0])),[3655.52723])
            axF.plot(f_window, self.baseLineRemove(Y_window[0]/np.amax(Y_window[0])), 'k', clip_on=True, label=label)
            self.result['amplitude'][gas] = self.baseLineRemove(Y_window[0]/np.amax(Y_window[0]))
            #axF.plot(f_window, Y_window_re[0]/(np.amax(Y_window_re[0])-np.amin(Y_window_re[0])), label=label+'_re')
            #axF.plot(f_window, Y_window_im[0]/(np.amax(Y_window_im[0])-np.amin(Y_window_im[0])), label=label+'_im')
            #axF.plot(f_window, self.baseLineRemove(Y_window_re[0]/np.amax(Y_window[0])))#, label=label+'_re')
            #axF.plot(f_window, self.baseLineRemove(Y_window_im[0]/np.amax(Y_window[0])))#, label=label+'_im')
            if gas == 'Ch2':
                pass
            else:
                pass
                axP.errorbar(f_window,P_window[0],yerr=P_window[1], color='r', ecolor='r',linewidth=2,label='Phase')
                #axP.axhline(y=0.5, color='b', linestyle='--', label='Phase_aid')
                #axP.axhline(y=1.5, color='b', linestyle='--')
            #plot(f_window,P_window,'r')
            #axF.set_ylim([0,0.2])
            #axF.set_xlim([200,4500])
            axP.set_yticks([-1,-0.5,0,0.5,1])
            axP.grid(visible=True,linestyle='--',linewidth='0.3',c='b')
            axP.set_ylim([-1.3,1.3])
            axP.legend(loc=(0.02,0.6),ncol=2,fontsize=10)
            #axF.set_yscale('log')
            #axF.set_ylim([10**-2.5,10**0])
            axF.legend(loc=(0.02,0.8),ncol=2,fontsize=10)
    
            i=i+1
            if ifsaveFFT:
                self.wksFFT.from_list(0, f_window, lname="Frequency", axis='X')
                self.wksFFT.from_list(i*3-2, self.baseLineRemove(Y_window[0]/np.amax(Y_window[0])), lname=label, axis='Y')
                #self.wks.from_list(i*3-1, P_window[0], lname=label, axis='Y')
                #self.wks.from_list(i*3, P_window[1], lname=label, axis='E')
            if ifsaveT:
                #self.wks.from_list(0, self.rebin_factor(self.delayB,10), 'X')
                #self.wks.from_list(i, self.rebin_factor((self.specBigBottleB[gas]-np.mean(self.specBigBottleB[gas])),10), lname=gas, axis='Y')
                self.wksT.from_list(0, self.smoothedT['T']*1e15, 'X')
                self.wksT.from_list(i, self.smoothedT[gas], lname=gas, axis='Y')

        fig.tight_layout()
        save_obj(self.result, pl.PureWindowsPath(self.savePath, self.folder+'_result'+r'.pkl'))
        save_obj(self.fftSB, pl.PureWindowsPath(self.savePath, self.folder+'_fftSB'+r'.pkl'))
        #plt.savefig(os.path.join(os.path.join(self.savePath,r'fft_logy.png')),dpi=720,bbox_inches='tight',pad_inches=0,transparent=True)
        #plt.savefig(os.path.join(os.path.join(self.savePath,r'fft.png')),dpi=720,bbox_inches='tight',pad_inches=0,transparent=True)
        plt.show()


    def show_Spectra(self, ifsaveT=False):
        gs = gridspec.GridSpec(2, 3)
        #gs = gridspec.GridSpec(1, 1)
        fig = plt.figure(figsize=(20,8))
        #ax = fig.add_subplot(111)
        lab=['Mass1','Mass2','Mass16','Mass17','Mass18','Ch2+Ch4','Mass1+Mass16', 'Mass1+Mass17', 'Mass16+Mass17','Mass17+Mass18']
        i=0
        if ifsaveT:
            op.set_show(show=True)
            self.wks = op.new_sheet('w',lname=str('Time')+str('_')+self.molecule+str('_')+self.folder)
        for gas in self.specBigBottleB.keys():
            print(gas)
            if 'filter' in gas or 'window' in gas or 'rebin' in gas or 'com' in gas:
                continue
            if gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
                pass
            else:
                print('Skip!')
                continue
            label=lab[i]
            
            ax = fig.add_subplot(gs[i])
            i=i+1
            delay = self.delayB
            rebindelay = self.rebin_delay
            #else:
            #    delay = self.rebin_delay
            ax.plot(delay*10**15,
                     self.specBigBottleB[gas]/(np.amax(self.specBigBottleB[gas][5500:6000])-np.amin(self.specBigBottleB[gas][5500:6000])), label=gas)#/np.amax(self.specBigBottleB[gas])
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
            if ifsaveT:
                self.wks.from_list(0, delay*10**15, 'X')
                self.wks.from_list(i, self.specBigBottleB['filter_'+gas], lname=label, axis='Y')
        #plt.legend()
        fig.tight_layout()
        plt.show()

    def calDrift(self, _cRange, gas='Ch8'):
        plt.clf()
        _iS = np.array(np.zeros(self.specBigBottleB[gas].shape))
        for i in range(self.specBigBottleB[gas].shape[0]):
            #_iS[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], (3643.52883-0.1)/33.35641*1e12, (3643.52883+0.1)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            _iS[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], 100/33.35641*1e12, 4000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
        _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=2000)
        indexMax = []
        for i in range(self.specBigBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange)
            plt.plot(_delayRange,_inter,label=str(i))
            #indexMax = indexMax + [sps.argrelextrema(_inter, np.greater)[0][-1]]
            indexMax = indexMax + [np.argmax(_inter)]
        plt.xlabel('Delay (fs)')
        plt.ylabel('a.u.')
        plt.legend()
        plt.show()

        _ref = indexMax[int(self.specBigBottleB[gas].shape[0]/2)]
        _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/2000
        for i in range(self.specBigBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange+_shift[i])
            plt.plot(_delayRange,_inter, label=str(i))
        plt.xlabel('Delay (fs)')
        plt.ylabel('a.u.')
        plt.legend()
        plt.show()
        return _shift
    
    def calDrift2(self, _cRange, gas='Ch8'):
        '''
        calibrate by drift of the strech mode oscillation
        '''
        plt.clf()
        _iS = np.array(np.zeros(self.specBigBottleB[gas].shape))
        #_iS2 = np.array(np.zeros(self.specBigBottleB[gas].shape))
        for i in range(self.specBigBottleB[gas].shape[0]):
            #_iS[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], (3655.52723-20)/33.35641*1e12, (3655.52723+20)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            _iS[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], 100/33.35641*1e12, 4000/33.35641*1e12, 1/self.stepSize)
            #_iS2[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], 100/33.35641*1e12, 4000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
        #iS2 = sp.interpolate.RectBivariateSpline(range(_iS2.shape[0]), self.delayB*1e15, _iS2)
        _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=20000)
        indexMax = []
        for i in range(self.specBigBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange)
            #plt.plot(_delayRange,_inter)
            #indexMax = indexMax + [sps.argrelextrema(_inter, np.greater)[0][-1]]
            indexMax = indexMax + [np.argmax(np.abs(_inter))]
            
        print(indexMax)
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        #plt.show()
        #_ref = sum(indexMax[int(self.specBigBottleB[gas].shape[0]/2)-5:int(self.specBigBottleB[gas].shape[0]/2)+5])/10
        _ref =indexMax[int(self.specBigBottleB[gas].shape[0]/2)]
        _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/20000
        for i in range(self.specBigBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange+_shift[i])
            #delayRangeA = np.linspace(start=50,stop=100, num=2000)#Used to compare the zreo delay when calibrating using stretch mode frequency to calibrate
            #_inter2 = iS2.ev(i,delayRangeA+_shift[i])#
            #plt.plot(_delayRange,_inter)
            #plt.plot(delayRangeA,_inter2)
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        #plt.show()
        return _shift

    #def delayCorrection(self, _cRange = [300,304.5]):def remo
    def delayCorrection(self, _cRange = [30,90]):
    #def delayCorrection(self, _cRange = [900,1300]):
        xxx = 0
        while True:
            try:
                #_shift = self.calDrift(_cRange = _cRange)
                _shift = self.calDrift2(_cRange = _cRange)
                break
            except IndexError:
                print('Wrong Region! Try again!')
                if xxx>10:
                    _cRange = [np.array(_cRange)[0], np.array(_cRange)[1]+0.5]
                _cRange = np.array(_cRange)+0.5
                xxx = xxx+1
        iinter = {}
        for gas in self.specBigBottleB.keys():
            iS = sp.interpolate.RectBivariateSpline(range(self.specBigBottleB[gas].shape[0]), self.delayB_noPad*1e15, self.specBigBottleB[gas])
            inter = np.zeros(self.specBigBottleB[gas].shape)
            for i in range(self.specBigBottleB[gas].shape[0]):
                inter[i] = iS.ev(i, self.delayB_noPad*1e15+_shift[i])
            iinter[gas] = inter
        self.specBigBottleB = iinter
    def wlt(self, gas):
        PHztocm = 33356.40952
        filt_freqs_in_cm = np.array([812,1344,2155,3657])
        bw = 60
        wavelet = 'cmor3-1'
        signal = self.specBigBottleB['window_'+gas] - np.mean(self.specBigBottleB['window_'+gas])
        signal = sps.savgol_filter(signal, window_length=50, polyorder=1,deriv=1,delta=20, mode='nearest')
        signal = signal/np.max(signal)
        delay = self.delayB*1e15
        delay = delay[np.where(delay>200)]
        scales = np.logspace(1.8, 3, num=100)
        dt = delay[1] - delay[0]

        frequencies = pywt.scale2frequency(
            wavelet, scales) / (delay[1] - delay[0])
        print('min freq (cm-1): ' + str(frequencies.min()*PHztocm))
        print('max freq (cm-1): ' + str(frequencies.max()*PHztocm))

        [coefficients, frequencies] = pywt.cwt(
            signal, scales, wavelet, delay[1] - delay[0])
        amplitude = np.abs(coefficients)
        amplitude = amplitude/np.max(amplitude)
        amplitude = amplitude[:,-np.size(delay):]
        period = 1/frequencies
        coi = np.sqrt(2)*6/(2*np.pi)/frequencies*2
        fig, ax = plt.subplots(
            2+filt_freqs_in_cm.size,
            1,
            gridspec_kw={'height_ratios': [3, 1, ] + filt_freqs_in_cm.size*[1]},
            sharex=True)
        fig.set_figheight(8)
        fig.set_figwidth(8)
        im = ax[0].pcolormesh(
            delay,
            frequencies*PHztocm,
            #10*np.log10(amplitude),
            amplitude,
            cmap='jet',
            vmax = 0.05)

        ax[0].plot(delay[-1]-coi,frequencies*PHztocm,color='white', linestyle='dashed',linewidth=3)
        ax[0].plot(coi,frequencies*PHztocm,color='white', linestyle='dashed',linewidth=3)

        ax[0].set_xlim([delay.min(), delay.max()])
        ax[0].set_ylim([0, 5000])
        ax[0].set_ylabel('frequency [cm-1]')
        ax[0].xaxis.set_ticklabels([])
        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes('top', size='7%', pad='2%')
        cb = fig.colorbar(im, cax=cax, orientation='horizontal')
        cb.set_label('amplitude')
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
            cc=int(2*np.sqrt(2)*6/(2*np.pi)/np.amin(frequencies[np.where(np.abs(filt_freq_in_cm-frequencies*PHztocm)<40)])/dt)
            occ = amplitude[filter, :].sum(axis=0)
            occ = occ[:-cc]/norm
            ax[1+idx].plot(
                delay[:-cc],
                #10*np.log10(occ),
                occ[-np.size(delay):],
                #label='ampl. [dB], freq = ' + str(filt_freq_in_cm))
                label='freq = ' + str(filt_freq_in_cm))
            ax[1+idx].set_xlim([delay.min(), delay.max()])
            #ax[1+idx].set_ylim([-30, 1])
            ax[1+idx].xaxis.set_ticklabels([])
            ax[1+idx].legend(frameon=0, loc='upper right')

        ax[-1].plot(delay, signal[-np.size(delay):])
        ax[-1].set_xlim([delay.min(), delay.max()])
        ax[-1].set_ylim([
            np.mean(signal)-2*np.std(signal),
            np.mean(signal)+2*np.std(signal)])
        ax[-1].set_ylabel('mass amplitude')
        ax[-1].set_xlabel('time (fs)')

        # fig.tight_layout()
        #plt.savefig('wavelet_transform_mass_' + str(gas) + '.png')
        plt.show()
    def cal_ratio(self):
        self.ratio = {}
        #if volts out=1V, volts in: 
        #1V, 10mV, 10mV, 1V ,5V and 10mV #4.5E+14_H2O 9e-7 mbar
        #1V, 20mV, 10mV, 1V ,5V and 20mV #7.2E+14_H2O 5e-7 mbar
        #1V, 20mV, 10mV, 2V ,5V and 20mV #8.9E+14_H2O 5e-7 mbar
        if self.folder == r'4.5E+14_H2O':
            boxcarA=[1,0.01,0.01,1,5,0.01]
            u = np.sum(self.specBigBottleB['Ch4'][-9000:])*0.01#*9/5 #9/5 is to compensate the pressure
        elif self.folder == r'7.2E+14_H2O':
            boxcarA=[1,0.02,0.01,1,5,0.02]
            u = np.sum(self.specBigBottleB['Ch4'][-9000:])*0.01
        elif self.folder == r'8.9E+14_H2O':
            boxcarA=[1,0.02,0.01,2,5,0.02]
            u = np.sum(self.specBigBottleB['Ch4'][-9000:])*0.01
        #print(u)
        i=0
        for gas in self.specBigBottleB.keys():
            self.ratio[gas] = np.sum(self.specBigBottleB[gas][-9000:])/u*boxcarA[i]
            if gas == 'Ch2':
                self.ratio[gas] = int(self.ratio[gas]*10)/10
            else:
                self.ratio[gas] = int(self.ratio[gas])
            i=i+1

if __name__ == '__main__':
    for ff in [r'3.6E+14_H2',r'6.0E+14_H2',r'8.0E+14_H2']:#r'3.6E+14_H2',r'6.0E+14_H2',r'8.0E+14_H2'
        d = FFT_ionS(ff)
        if d.checkSavedData():
            d.read()
            d.delayCorrection()
        d.transition()
        #d.cal_ratio()
        #print(d.ratio)
        d.findZeroDelay3()
        d.show_Spectra(ifsaveT=False)
        #d.FFT3(windowSize=100, delayRange=[300*1E-15,1000*1E-15], rebinF=1,paddingF = 5, useWindow=True, zeroDirection='left', phaseCompensate=False, smooth=True,test = False)
        d.FFT3(windowSize=90, delayRange=False, rebinF=1,paddingF = 10, useWindow=True, zeroDirection='left', phaseCompensate=True, smooth=True,test = False)
        d.show_FFT(ifsaveFFT=False,ifsaveT=False)