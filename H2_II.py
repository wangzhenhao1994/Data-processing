#calibrate the delay drift for different measurement, cut the data if necessray
#do FFT and show FFT and phase
#Main function if self.show_FFT()
#smooth function doesn't influence the phase, checked

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
from skimage.restoration import denoise_wavelet
import math

my_path = os.path.abspath(__file__)

mpl.rcParams['lines.linewidth'] = 1
plt.rcParams["font.family"] = "arial"
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.it'] = 'Arial:italic'
mpl.rcParams['mathtext.rm'] = 'Arial'


SMALL_SIZE = 16
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
    def __init__(self,folder,caliFolder=''):
        '''

        '''
        self.folder = folder
        self.caliFolder = caliFolder
        self.filename = filename
        self.saveRef = str(filename)
        self.filepath = []
        self.delayB,self.stepSize = np.linspace(start=0,stop=1300, num=13000,endpoint=False,retstep=True)
        self.delayB = self.delayB*10**-15
        self.stepSize = self.stepSize*10**-15
        self.delayB_noPad = np.linspace(start=0,stop=1299.9, num=13000)*10**-15
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
        self.fftSB_H2 = {}
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

        #setting for data from digitizer
        #self.channelSize = 12032#24000#1536
        #self.scanLengthD = 3320#1200
        #self.peakRange = [-100, 100]  # range of the peak
        #self.delayD = np.arange(self.scanLengthD)/self.scanLengthD*100*2*2*3.33564*10**-15#*3602/3647
        #self.spectraBottleD = {}
        #self.fftSD = {}
        #self.stftSD = {}
        #self.dataD = 0

        self.rootPath= pl.PureWindowsPath(r'D:\SF_FFT')
        self.savePath= pl.PureWindowsPath(os.path.join(r'D:\SF_FFT\processedData\inter',self.molecule,folder))

    def read(self):
        totalCounter = {}
        #for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
        #    self.fftSB[gas] = 0
        for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
            self.specBigBottleB[gas] = np.zeros((1000,13000))
            totalCounter[gas] = 0
        #if not os.path.exists(os.path.join(self.savePath,r'plotData')):
        #    os.mkdir(os.path.join(self.savePath,r'plotData'))
        for file in os.listdir(self.savePath):
            filename = os.fsdecode(file)
            print(filename)
            if filename.endswith(".pkl"):
                self.interSpectraBottleB = load_obj(os.path.abspath(os.path.join(self.savePath, filename)))
                for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
                    print(self.interSpectraBottleB[gas].shape)
                    for i in range(self.interSpectraBottleB[gas].shape[0]):
                        self.specBigBottleB[gas][totalCounter[gas]] = self.interSpectraBottleB[gas][i]
                        totalCounter[gas] = totalCounter[gas]+1
        print(totalCounter)

        for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
            self.specBigBottleB[gas] = self.specBigBottleB[gas][:totalCounter[gas]]
        zeroIndex2 = self.findZeroDelay2()
        zeroIndex2Min = np.amin(zeroIndex2)
        #zeroIndex2Max = np.amax(zeroIndex2)
        for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
            for i in range(self.specBigBottleB[gas].shape[0]):
                self.specBigBottleB[gas][i] = np.concatenate((self.specBigBottleB[gas][i][(zeroIndex2[i]-zeroIndex2Min):], np.full(zeroIndex2[i]-zeroIndex2Min,self.specBigBottleB[gas][i][-1])))
        print(self.specBigBottleB['Ch0'].shape)
    def interFFT(self, y):
        n = len(y)
        delta = self.stepSize#/(self.interNum+1)
        self.dw = 1/((n)*self.stepSize)/1E12*33.35641
        #print('self.dw is '+str(self.dw))
        f = np.fft.rfftfreq(n, delta)/1E12*33.35641  # frequency unit cm-1
        f = f[f>0]
        #y=np.fft.fftshift(y)
        fft_y = np.fft.rfft(y)

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
        self.ifDelayRange = delayRange
        self.delayRange4FFT = delayRange
        for gas in self.phaseSpecBottleB.keys():
            self.smoothedT[gas]=0
            for i in range(self.phaseSpecBottleB[gas].shape[0]):
                #becareful, here the sign of the amplitude is changed
                interSpec = -self.phaseSpecBottleB[gas][i][-self.specBigBottleB[gas].size:]
                #becareful, here the sign of the amplitude is changed
                _,interSmoothedT,self.smoothedT['T'] = self.inter_window(interSpec,self.delayB,windowSize=0,direction=zeroDirection, useWindow=False, phaseCompensate=False,smooth = False)
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
                #plt.plot(self.refSpectra4AbsPhase)
                #plt.plot(y)
                #plt.show()
                if rebinF < 1.5:
                    pass
                else:
                    y = self.rebin_factor(y,rebinF)
                    t = self.rebin_factor(t,rebinF)
                #self.interNum = 100
                #y = self.interInterp(t,y,self.interNum)
                self.inputLength = np.size(y)
                #print('Input length is '+str(self.inputLength)+' !')
                f, Y = self.interFFT(y)
                self.L=np.size(f)
                if i == 0:
                    _interY = np.zeros((self.phaseSpecBottleB[gas].shape[0],Y.size),dtype=np.complex64)
                _interY[i] = Y
            self.fftSB['window_'+gas+'_fft'] = _interY
            self.fftSB[gas] = 0
            __, self.fftSB['ref4Phase'] = self.interFFT(self.refSpectra4AbsPhase)

        self.fftSB['window_'+'fre'] = f
        #print(self.folder)
        self.smooth = smooth
        save_obj(self.fftSB,os.path.join(self.savePath,str(smooth)+str(self.folder)+r'_fftSB.pkl'))

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
            n=5
            print(self.specBigBottleB[gas].shape[0],n)
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
        for gas in list(self.spectraBottleD):
            if 'rebin' not in gas:
                self.spectraBottleD[gas], self.spectraBottleD['window_'+gas], self.delayD = self.inter_window(self.spectraBottleD[gas], self.delayD, windowSize=windowSize, direction=direction, useWindow=useWindow)
            else:
                self.spectraBottleD[gas], self.spectraBottleD['window_rebin_'+gas], self.rebin_delay = self.inter_window(self.spectraBottleD[gas], self.rebin_delay, windowSize=windowSize, direction=direction, useWindow=useWindow)
                self.windowSize = windowSize
    
    def window2(self, windowSize=0, direction='left'):
        '''
        windowSize is in fs.
        '''
        for gas in list(self.specBigBottleB):
            if 'rebin' not in gas:
                self.specBigBottleB[gas], self.specBigBottleB['window_'+gas], self.delayB = self.inter_window2(self.specBigBottleB[gas], self.delayB, windowSize=windowSize, direction=direction)
            else:
                self.specBigBottleB[gas], self.specBigBottleB['window_rebin_'+gas], self.rebin_delay = self.inter_window2(self.specBigBottleB[gas], self.rebin_delay, windowSize=windowSize, direction=direction)
        for gas in list(self.spectraBottleD):
            if 'rebin' not in gas:
                self.spectraBottleD[gas], self.spectraBottleD['window_'+gas], self.delayD = self.inter_window2(self.spectraBottleD[gas], self.delayD, windowSize=windowSize, direction=direction)
            else:
                self.spectraBottleD[gas], self.spectraBottleD['window_rebin_'+gas], self.rebin_delay = self.inter_window2(self.spectraBottleD[gas], self.rebin_delay, windowSize=windowSize, direction=direction)
    
    def inter_window(self, data, delay, windowSize=0, direction='left', useWindow = True, phaseCompensate = True, smooth = True):
        '''
        windowSize is in fs.
        '''

        if self.ifDelayRange:#
            print('Using delay range!')
            data = data[:int(np.size(data)/2)*2]
            __len = np.size(data)
            delay = delay[:__len]
            windowSize = int(self.delayRange4FFT[0]/self.stepSize)
            #print('window size is '+ str(windowSize)+str(' !'))
            self.windowSize = windowSize
            if smooth:
                data=denoise_wavelet(data, sigma=5, wavelet='sym5', wavelet_levels=None)
                #data=sps.savgol_filter(data, window_length=10, polyorder=1,deriv=1,delta=10, mode='nearest')
            window=polynomial(data, order=1, plot=False)
            refDelay = np.arange(np.size(window))*self.stepSize*1e15
            omega=[588.26682,1346.61947,1479.61893,1588.32041,1704.69493,1821.06945,1945.11702,2053.8185,2189.37563,3294.29418,3530.87974,3771.30183,4161.34831,357.69977,815.04447,4518.73414,4718.93935,4911.44435,5110.36619,5463.29203,3081.36344]
            oscA=np.sum(np.array([np.sin(2*np.pi*refDelay/33356.40952*f) for f in omega]),0)
            self.refSpectra4AbsPhase = np.hanning(np.size(window))*oscA[-(__len-windowSize):]
            if useWindow:
                #window2=self.apply_triwindow(window)
                window2=self.apply_hannwindow(window)
                #window2=self.apply_hammwindow(window)
                #plt.plot(self.refSpectra4AbsPhase)
                #plt.show()
            else:
                window2=window


            if phaseCompensate:
                window2 = np.append(np.zeros(windowSize),window2,axis=0)
                window = np.append(np.zeros(windowSize),window,axis=0)
                self.refSpectra4AbsPhase = np.append(np.zeros(windowSize),self.refSpectra4AbsPhase,axis=0)
        else:
            data = data[:int(np.size(data)/2)*2]
            
            __len = np.size(data)
            delay = delay[:__len]
            windowSize = int(windowSize*1e-15/self.stepSize)
            #print('window size is '+ str(windowSize)+str(' !'))
            self.windowSize = windowSize
            refDelay = delay*1e15#np.arange(np.size(window))*self.stepSize*1e15
            omega=[588.26682,1346.61947,1479.61893,1588.32041,1704.69493,1821.06945,1945.11702,2053.8185,2189.37563,3294.29418,3530.87974,3771.30183,4161.34831,357.69977,815.04447,4518.73414,4718.93935,4911.44435,5110.36619,5463.29203,3081.36344]
            oscA=np.sum(np.array([np.sin(2*np.pi*refDelay/33356.40952*f) for f in omega]),0)
            #data=self.butter_bandpass_filter(data, 300/33.35641*1e12, 18000/33.35641*1e12, 1/self.stepSize)
            self.refSpectra4AbsPhase =oscA#self.butter_bandpass_filter(oscA, 300/33.35641*1e12, 18000/33.35641*1e12, 1/self.stepSize) 
            if smooth:
                #pass
                #data=denoise_wavelet(data, sigma=25, wavelet='sym5', wavelet_levels=None)
                a,b,c,d = [20,1,1,3]
                data=sps.savgol_filter(data, window_length=a, polyorder=b,deriv=c,delta=d, mode='nearest')
                self.refSpectra4AbsPhase=sps.savgol_filter(self.refSpectra4AbsPhase, window_length=a, polyorder=b,deriv=c,delta=d, mode='nearest')

            if direction == 'left':
                window = data[-(__len-windowSize):]
                self.refSpectra4AbsPhase = self.refSpectra4AbsPhase[-(__len-windowSize):]
                if phaseCompensate:
                    pass
                else:
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
            
            window = polynomial(window, order=15, plot=False)

            #window =self.baseLineRemove(window)#shift the base line to zero


            #if not math.log2(np.size(window)).is_integer():
            #    #print('Fill the data to length if power of 2!')
            #    window = np.append(window,np.zeros(2**(math.ceil(math.log2(np.size(window))))-np.size(window)),axis=0)
            if useWindow:
                #window2=self.apply_triwindow(window)
                window2=self.apply_hannwindow(window)
                self.refSpectra4AbsPhase = self.apply_hannwindow(self.refSpectra4AbsPhase)
                #window2=self.apply_hammwindow(window)
                #plt.plot(self.refSpectra4AbsPhase)
                #plt.show()
            else:
                window2=window


            if phaseCompensate:
                window2 = np.append(np.zeros(windowSize),window2,axis=0)
                window = np.append(np.zeros(windowSize),window,axis=0)
                self.refSpectra4AbsPhase = np.append(np.zeros(windowSize),self.refSpectra4AbsPhase,axis=0)
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
            #np.arange(delay[0]-(paddingSize+1) *
            #          delayStep, delay[0], delayStep),
            delay,
            np.arange(paddingSize)*delayStep+delay[-1]
        ))
        data = np.concatenate(
            #(np.zeros(paddingSize),inter_data, np.zeros(paddingSize)), axis=0)
            (inter_data, np.zeros(paddingSize)), axis=0)
        self.refSpectra4AbsPhase = np.concatenate(
            #(np.zeros(paddingSize),inter_data, np.zeros(paddingSize)), axis=0)
            (self.refSpectra4AbsPhase, np.zeros(paddingSize)), axis=0)
        data = data[:int(np.size(data)/2)*2]
        delay = delay[:len(data)]
        return data, delay

    def interpS(self, interNum=1):
        interpDelay = np.arange(
            self.delayB[0], self.delayB[-1], (self.delayB[1]-self.delayB[0])/(interNum+1))
        for gas in self.specBigBottleB.keys():
            iS = np.interp(interpDelay, self.delayB, self.specBigBottleB[gas])
            self.specBigBottleB[gas] = iS
            # #print(np.shape(iS))
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

    #def runingAverage(self, n=5):
    #    def runAve(x, N): return np.convolve(x, np.ones(N)/N, mode='valid')
    #    for gas in self.specBigBottleB.keys():
    #        self.specBigBottleB[gas] = runAve(self.specBigBottleB[gas], n)
    #        new_size = len(self.specBigBottleB[gas])
    #    self.delayB = self.delayB[:new_size]
#
    #def smooth(self, windowSize=100, order=3):
    #    for gas in self.specBigBottleB.keys():    
    #        self.specBigBottleB[gas]=sps.savgol_filter(self.specBigBottleB[gas], windowSize, order) # window size 51, polynomial order 3
    #
    #def smooth2(self,data,windowSize=100,order=3):
    #    return sps.savgol_filter(data, windowSize, order)

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
            _Spec = self.butter_bandpass_filter(_Spec, 100/33.35641*1e12, 1000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            #plt.plot(_Spec)
            zeroIndex =zeroIndex + [(_Spec[:2000]).argmax()]
        #plt.show()
        return zeroIndex

    def findZeroDelay2(self):
        zeroIndex = []
        for i in range(self.specBigBottleB['Ch2'].shape[0]):
            _Spec = np.abs(self.specBigBottleB['Ch2'][i])+np.abs(self.specBigBottleB['Ch0'][i])
            _Spec = self.butter_bandpass_filter(_Spec, 100/33.35641*1e12, 1000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            plt.plot(_Spec)
            
            zeroIndex =zeroIndex + [np.abs((_Spec[:2000])).argmax()]
        plt.show()
        return zeroIndex

    def findZeroDelay3(self):
        zeroIndex = []
        _Spec = np.abs(self.specBigBottleB['Ch0'])/np.max(np.abs(self.specBigBottleB['Ch0']))#+np.abs(self.specBigBottleB['Ch8'])/np.max(np.abs(self.specBigBottleB['Ch8']))
        _Spec = self.butter_bandpass_filter(_Spec, 300/33.35641*1e12, 4500/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        zeroIndex =np.abs((_Spec[:2000])).argmax()
        self.zeroIndex = zeroIndex # this number will be used for phase correction
        __inter = self.specBigBottleB_noPad.copy()
        self.specBigBottleB_noPad = {}
        for gas in self.specBigBottleB.keys():
            self.specBigBottleB[gas] = self.specBigBottleB[gas][zeroIndex:]
            self.specBigBottleB_noPad[gas]={}
            for j in range(__inter[gas].shape[0]):
                self.specBigBottleB_noPad[gas][j]=-1*__inter[gas][j][zeroIndex:] #change the sign of the amplitude
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
        self.interfftSB = load_obj(os.path.abspath(os.path.join(d.savePath, str(self.smooth)+str(self.folder)+r'_fftSB.pkl')))
        self.dcRange=0
        fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True,figsize=(10,10))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel("Frequency ($\mathrm{cm^{-1}}$)")
        plt.ylabel("Amplitude (a.u.)", labelpad=25)
        plt.text(1.18,0.35,'Relative Phase v.s. Mass 2 ($\pi$)', rotation = 270)
        axPP = []
        for axx in ax.flatten():
            axPP = axPP + [axx.twinx()]
        axPP[0].get_shared_y_axes().join(axPP[0], axPP[1])

        lab=['Mass 2','Mass 1']
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
        boxcarA['Ch0'],boxcarA['Ch2'] = [1,0.01]

        for gas in ['Ch2','Ch0']:
            if 'filter' in gas  or 'window' in gas or 'rebin' in gas or 'com' in gas:
                continue
            elif gas not in ['Ch2','Ch0']:
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
            #if self.smooth and self.folder==r'7.2E+14_H2O':
            #    P_inter = np.angle(self.interfftSB[_preP+gas+'_fft'])
            #else:
            P_inter = np.angle(self.fftSB[_preP+gas+'_fft'])
            
            if gas == 'Ch2':
                P_ref = P_inter
                P_window =  np.array([np.mean(P_inter,axis=0)-np.angle(self.fftSB['ref4Phase']), np.std(P_inter,axis=0)])#np.mean(P_inter,axis=0)-np.angle(self.fftSB['ref4Phase'])
                #plt.plot(f_window,np.angle(self.fftSB['ref4Phase']))
                #plt.plot(f_window,np.unwrap(np.mean(P_inter,axis=0)))
                #plt.plot(np.arctan2(Y_window_im[0],Y_window_re[0]))
                #plt.xlim([3500,3600])
                #plt.show()
            else:
                P_inter = P_inter-P_ref
                for iii in range(10):
                    P_inter.sort(axis=0)
                    P_inter2=P_inter.copy()
                    P_inter3=P_inter.copy()
                    P_inter2[-1]=P_inter2[-1]-2*np.pi
                    P_inter3[0]=P_inter3[0]+2*np.pi
                    P_inter[-1]=np.where(np.std(P_inter,axis=0)<np.std(P_inter2,axis=0),P_inter[-1],P_inter2[-1])
                    P_inter[0]=np.where(np.std(P_inter,axis=0)<np.std(P_inter3,axis=0),P_inter[0],P_inter3[0])

                #P_window =  np.array([np.mean(P_inter,axis=0), np.sqrt(np.std(P_inter,axis=0)**2+np.std(P_ref,axis=0)**2)])
                P_window =  np.array([np.mean(P_inter,axis=0), np.sqrt(np.std(P_inter,axis=0)**2+0)])

            
            aa = len(f_window[(f_window<100)])
            bb=len(f_window[(f_window<20000)])
            #aa = len(f_window[(f_window<100)])
            #bb=len(f_window[(f_window<50000)])
            f_window = f_window[aa:bb]
            #print('The step of frequency is '+str(np.mean(f_window[1:10]-f_window[0:9])))
            Y = Y[:,aa:bb]
            Y_window = Y_window[:,aa:bb]
            Y_window_im = Y_window_im[:,aa:bb]
            Y_window_re = Y_window_re[:,aa:bb]
            P_window = P_window[:,aa:bb]
            P_window[0]=np.where(P_window[0]<-np.pi+1,P_window[0]+2*np.pi,P_window[0])
            P_window[0]=np.where(P_window[0]>np.pi-1,P_window[0]-2*np.pi,P_window[0])
            P_window = P_window/np.pi
            self.result['frequency'] = f_window
            self.result['phase'][gas] = P_window
            omega=[588.26682,1346.61947,1479.61893,1588.32041,1704.69493,1821.06945,1945.11702,2053.8185,2189.37563,3294.29418,3530.87974,3771.30183,4161.34831,357.69977,815.04447,4518.73414,4718.93935,4911.44435,5110.36619,5463.29203,3081.36344,5684]

                
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
                inter = inter + np.where(np.abs(f_window-w)<20,P_window,0)
            inter = np.where(inter==0,np.inf,inter)
            P_window = inter
            
            _interF = self.baseLineRemove(sps.savgol_filter(Y_window[0], window_length=20, polyorder=0,deriv=0,delta=2, mode='nearest'))
            axF.plot(f_window, _interF/np.max(_interF), 'k', clip_on=True, label=label)
            #self.result['amplitude'][gas] = self.baseLineRemove(Y_window[0]/np.amax(Y_window[0]))
            #axF.plot(f_window, Y_window_re[0]/(np.amax(Y_window_re[0])-np.amin(Y_window_re[0])), label=label+'_re')
            #axF.plot(f_window, Y_window_im[0]/(np.amax(Y_window_im[0])-np.amin(Y_window_im[0])), label=label+'_im')
            #axF.plot(f_window, self.baseLineRemove(Y_window_re[0]/np.amax(Y_window[0])))#, label=label+'_re')
            #axF.plot(f_window, self.baseLineRemove(Y_window_im[0]/np.amax(Y_window[0])))#, label=label+'_im')
            axP.errorbar(f_window,P_window[0],yerr=P_window[1], color='r', ecolor='r',linewidth=2,label='Phase')
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
                self.wksFFT.from_list(i*3-2, _interF/np.max(_interF), lname=label, axis='Y')
                self.wksFFT.from_list(i*3-1, P_window[0], lname=label, axis='Y')
                self.wksFFT.from_list(i*3, P_window[1], lname=label, axis='E')
                #self.wksFFT.from_list(i*3-3, f_window, lname="Frequency", axis='X', units = 'wavenumber')
                #self.wksFFT.from_list(i*3-2, self.baseLineRemove(Y_window[0]/np.amax(Y_window[0])), lname=label, axis='Y')
                #self.wksFFT.from_list(i*3-2, Y_window_re[0], lname=label, axis='Y', comments = 're')
                #self.wksFFT.from_list(i*3-1, Y_window_im[0], lname=label, axis='Y', comments = 'im')

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

    #def show_FFTD(self):
    #    #plt.figure()
    #    #unused_f, ref = self.fftS['CH3OH+']
    #    i=12
    #    for gas in self.fftSD.keys():
    #        #if 'filter' not in gas:
    #        #    continue
    #        [f, Y, _] = self.fftSD[gas]
    #        #Y=np.where(np.logical_and(f>=10, f<=4500),Y,0)
    #        if gas == 'window_H2+':
    #            plt.plot(f[self.dcRange:], Y[self.dcRange:]/1000, label=r'$\mathrm{H_2^+}$')
    #        #if gas == 'H+':
    #        #    plt.plot(f[self.dcRange:], Y[self.dcRange:]/1000, label=r'$\mathrm{H^+}$')
    #        if gas == 'window_O+':
    #            plt.plot(f[self.dcRange:], Y[self.dcRange:]/1000, label=r'$\mathrm{O^+}$')
    #        if gas == 'window_O+H2':
    #            plt.plot(f[self.dcRange:], Y[self.dcRange:]/1000, label=r'$\mathrm{H_2^++O^+}$')
    #            #print("The max of the FFT is ",np.amax( Y[self.dcRange:]/1000))
    #        if gas == 'window_O-H2':
    #            plt.plot(f[self.dcRange:], Y[self.dcRange:]/1000, label=r'$\mathrm{H_2^+-O^+}$')
    #        else:
    #            #plt.plot(f[self.dcRange:], Y[self.dcRange:]/1000, label=gas)
    #            pass
    #        plt.xlabel(r'$\mathrm{Frequency (cm^{-1})}$')
    #        plt.ylabel('a.u.')
    #        plt.legend(loc=2)
    #        #plt.xlabel('Frequency/cm-1')
    #        #plt.ylim([0,np.max(Y[self.dcRange:1000])*3/2])
    #        plt.xlim([0,4500])
    #        i=i+1
    #        if ifsave:
    #            self.wks.from_list(7, f[self.dcRange:], 'X')
    #            self.wks.from_list(i, np.abs(Y[self.dcRange:]), lname=gas, axis='Y')
    #    #plt.legend()
    #    plt.tight_layout()
    #    d='_100_120'

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
            #print(gas)
            if 'filter' in gas or 'window' in gas or 'rebin' in gas or 'com' in gas:
                continue
            if gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
                pass
            else:
                #print('Skip!')
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
        standard = 0
        standard = self.specBigBottleB['Ch8']-self.specBigBottleB['Ch0']#-self.specBigBottleB['Ch6']
        _iS = np.array(np.zeros(self.specBigBottleB[gas].shape))
        #_iS2 = np.array(np.zeros(self.specBigBottleB[gas].shape))
        for i in range(self.specBigBottleB[gas].shape[0]):
            #_iS[i]=self.butter_bandpass_filter(standard[i], (3655.52723-20)/33.35641*1e12, (3655.52723+20)/33.35641*1e12, 1/self.stepSize)
            #_iS[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], 500/33.35641*1e12, 3800/33.35641*1e12, 1/self.stepSize)
            _iS[i]=self.butter_bandpass_filter(standard[i], 100/33.35641*1e12, 1500/33.35641*1e12, 1/self.stepSize)
            #_iS2[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], 100/33.35641*1e12, 4000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
        #iS2 = sp.interpolate.RectBivariateSpline(range(_iS2.shape[0]), self.delayB*1e15, _iS2)
        _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=20000)
        indexMax = []
        for i in range(self.specBigBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange)
            plt.plot(_delayRange,_inter)
            #indexMax = indexMax + [sps.argrelextrema(_inter, np.greater)[0][-1]]
            indexMax = indexMax + [np.argmax(np.abs(_inter))]
            
        #print(indexMax)
        plt.xlabel('Delay (fs)')
        plt.ylabel('a.u.')
        plt.show()
        #_ref = sum(indexMax[int(self.specBigBottleB[gas].shape[0]/2)-5:int(self.specBigBottleB[gas].shape[0]/2)+5])/10
        _ref =indexMax[int(self.specBigBottleB[gas].shape[0]/2)]
        _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/20000
        for i in range(self.specBigBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange+_shift[i])
            #delayRangeA = np.linspace(start=50,stop=100, num=2000)#Used to compare the zreo delay when calibrating using stretch mode frequency to calibrate
            #_inter2 = iS2.ev(i,delayRangeA+_shift[i])#
            plt.plot(_delayRange,_inter)
            #plt.plot(delayRangeA,_inter2)
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        plt.show()
        return _shift

    #def delayCorrection(self, _cRange = [300,307]):
    def delayCorrection(self, _cRange = [1,150]):
    #def delayCorrection(self, _cRange = [900,1300]):
        xxx = 0
        while True:
            try:
                #_shift = self.calDrift(_cRange = _cRange)
                _shift = self.calDrift2(_cRange = _cRange)
                break
            except IndexError:
                #print('Wrong Region! Try again!')
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

    def removeHydrogenS(self):
        savePathH2= pl.PureWindowsPath(os.path.join(r'C:\Users\user\Desktop\Data_newTOF\dataProcessing\H2',self.caliFolder))
        if os.path.exists(os.path.join(savePathH2,r'totalSpec.pkl')):
            self.specBigBottleB_H2 = load_obj(os.path.join(savePathH2,r'totalSpec.pkl'))
        for gas in self.specBigBottleB_H2.keys():
            self.specBigBottleB_H2[gas] = self.specBigBottleB_H2[gas].sum(axis = 0)

    def phaseRetrive(self,omega=[3655.52723]):
        #[526.49165,625.20883,699.99458,810.67748,882.47179,1145.71762,1343.15199,1594.43209,1696.1407,2153.82946,2324.34096,2677.32968,3212.79562,3655.52723]
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
                    #print(zeroIndex)
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
            if gas == 'Ch8':
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

    

    def phaseRetrive2(self,omega=[526.49165,625.20883,699.99458,810.67748,882.47179,1145.71762,1343.15199,1594.43209,1696.1407,2153.82946,2324.34096,2677.32968,3212.79562,3655.52723]):
        """
        Retrieve the absolute phase
        """
        
        self.interSpectraBottleB = {}
        self.phaseBottle = {}

        for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8']:
            _shape = self.specBigBottleB_noPad[gas].shape[0]
            _shape = _shape-_shape%4
            self.interSpectraBottleB[gas] = self.specBigBottleB_noPad[gas][:_shape]
            #self.interSpectraBottleB[gas] = self.interSpectraBottleB[gas].reshape(int(_shape/6), 6, self.interSpectraBottleB[gas].shape[1]).sum(axis=1)
            self.interSpectraBottleB[gas] = self.interSpectraBottleB[gas].reshape(int(_shape/4), 4, self.interSpectraBottleB[gas].shape[1]).sum(axis=1)
            self.phaseBottle[gas] = {}

        for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8']:
            _zeroIndex = self.findZeroDelay()
            #print(_zeroIndex)
            for j in range(self.interSpectraBottleB[gas].shape[0]):
                self.phaseBottle[gas][str(j)]={}
                interPhase=[]
                
                _,_ss,_d = self.inter_window(self.interSpectraBottleB['Ch8'][j],self.delayB_noPad,150)
                y,t = self.inter_padding( _ss,_d,paddingSize=100000)
                [f_window, Y_window, _] = self.interFFT(t, y)
                f_window = f_window
                Y_window = Y_window[:len(f_window)]
                amp = np.array([np.mean(np.where(np.logical_and(f_window>=ff-20, f_window<=ff+20),Y_window,0)) for ff in omega])
                amp = np.reshape(amp,(1,len(amp)))
                _,_s,_d = self.inter_window(self.interSpectraBottleB[gas][j],self.delayB_noPad,150)
                y,t = self.inter_padding( _s,_d,paddingSize=100000)
                [f_window, Y_window, _] = self.interFFT(t, y)  
                f_window = f_window
                osc = _s
                hwindow = np.hanning(len(_d))
                cal_len = int(33356.40952/520.50412/self.stepSize/1e15*1)#choose how much steps to calculate
                for i in range(cal_len):
                    step=self.delayB[1]-self.delayB[0]
                    oscA=np.sum(np.dot(amp,np.array([np.sin(2*np.pi*f*(_d-(_zeroIndex[j]+i)*step)) for f in np.array(omega)/33.35641*1e12])),0)
                    #sps.windows.flattop(len(window),sym=FALSE)
                    #plt.plot(self.delayB_noPad,hwindow*(oscA/100))
                    #plt.plot(self.delayB_noPad,osc)
                    #plt.show()
                    y,t=self.inter_padding(hwindow*(oscA/np.amax(oscA))*0.01+osc/np.amax(polynomial(osc, order=15, plot=False)), _d, paddingSize = 100000)
                    f_inter,Y,_ = self.interFFT(t,y)
                    #plt.plot(f_inter,Y)
                    #plt.xlim([300,4500])
                    #plt.show()
                    interPhase = interPhase+[Y]
                #plt.contourf(interPhase)
                #plt.show()
                f_inter = f_inter
                index0 = np.argmin(np.abs(f_inter-300))
                index1 = np.argmin(np.abs(f_inter-4500))
                if gas=='Ch8' and j==0:
                    self.interPhase = np.array(interPhase)[:,np.argmin(np.abs(f_inter-10)):np.argmin(np.abs(f_inter-15000))]
                    save_obj(self.interPhase, os.path.join(self.savePath, r'interPhase.pkl'))
#                #if ifsavePhase:
                #    self.wks = op.new_sheet('w',lname=str('Phase')+str('_')+str(gas)+str('_')+str('%.1E' % Decimal(self.intensity))+str('_')+self.saveRef)
                for f in omega:
                    self.phaseBottle[gas][str(j)][str(f)] = np.sum(np.array(interPhase)[:,np.argmin(np.abs(f_inter-f+5)):np.argmin(np.abs(f_inter-f-5))],1)
                    #tt = self.delayB[:cal_len]-_zeroIndex[j]*self.delayStep
                    #res = self.fit_sin(tt,self.phaseBottle[gas][str(j)][str(f)],f/33.35641*1e12)
                    #self.phaseBottle[gas][str(j)][str(f)+'_fit'] = [tt, res]
                    #A, w, p, c = self.phaseBottle[gas][str(j)][str(f)+'_fit'][1]
                    #fitfunc = lambda t: A * np.sin(w*t + p) + c
                    #if f == 3643.52883:
                    #    plt.plot(tt, fitfunc(tt), "r-", label="y fit curve", linewidth=2)
                    #    plt.plot(tt,self.phaseBottle[gas][str(j)][str(f)])
            ##plt.show()
        save_obj(self.phaseBottle, os.path.join(self.savePath, r'phaseBottle.pkl'))

    def phaseRetrive3(self):
        """
        Retrieve the relative phase to the phase in 'Ch8'
        """
        omega=[526.49165,625.20883,699.99458,810.67748,882.47179,1145.71762,1343.15199,1594.43209,1696.1407,2153.82946,2324.34096,2677.32968,3212.79562,3655.52723]
        self.interSpectraBottleB = {}
        self.phaseBottle = {}

        for gas in ['Ch8','Ch0','Ch2','Ch4','Ch6']:
            _shape = self.specBigBottleB_noPad[gas].shape[0]
            _shape = _shape-_shape%6
            self.interSpectraBottleB[gas] = self.specBigBottleB_noPad[gas][:_shape]
            self.interSpectraBottleB[gas] = self.interSpectraBottleB[gas].reshape(int(_shape/6), 6, self.interSpectraBottleB[gas].shape[1]).sum(axis=1)
            self.phaseBottle[gas] = {}
            for j in range(1,self.interSpectraBottleB[gas].shape[0]):
                #print(gas)
                self.phaseBottle[gas][str(j)]={}
                interPhase=[]
                
                _,_ss,_d = self.inter_window(self.interSpectraBottleB['Ch8'][j], self.delayB_noPad,100, useWindow = False)
                _,_s,_d = self.inter_window(self.interSpectraBottleB[gas][j], self.delayB_noPad,100, useWindow = False)

                cal_len = int(33356.40952/526.49165/self.stepSize/1e15)#choose how much steps to calculate
                oscAmp = polynomial(_s,order=5)
                osc = _s[:(len(_s)-cal_len)]/(np.amax(oscAmp[-500:])-np.amin(oscAmp[-500:]))
                hwindow = np.hanning(len(osc))
                
                for i in range(cal_len):
                    oscaAmp = polynomial(_ss,order=5)
                    oscA=_ss[i:(len(_ss)-cal_len+i)]/(np.amax(oscaAmp[-500:])-np.amin(oscaAmp[-500:]))
                    y,t=self.inter_padding(hwindow*(osc+oscA), _d, paddingSize = np.size(osc)*2)
                    f_inter,Y = self.interFFT(t,y)
                    #plt.plot(f_inter,Y)
                    #plt.xlim([300,4500])
                    #plt.show()
                    interPhase = interPhase+[Y]
                f_inter = f_inter
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
        save_obj(self.phaseBottle, os.path.join(self.savePath, r'relativePhase', r'phaseBottle.pkl'))

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
        #print('min freq (cm-1): ' + str(frequencies.min()*PHztocm))
        #print('max freq (cm-1): ' + str(frequencies.max()*PHztocm))

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

    def STFTS(self, gas, windowsize=0, ratio=1):
        '''
        windowsize is window size in fs 
        '''
        windowsize = int(windowsize/self.stepSize/1E15)
        fig = plt.figure(figsize=(8, 8))
        # Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
        # the size of the marginal axes and the main axes in both directions.
        # Also adjust the subplot parameters for a square plot.
        gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7),
                              left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.05, hspace=0.05)

        ax = fig.add_subplot(gs[1, 0])
        ax_histx = fig.add_subplot(gs[0, 0])
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

        #if 'rebin' not in gas or 'filter' not in gas:*
        if 'rebin' in gas:
            f, t, self.stftSB[gas] = sps.stft(self.specBigBottleB[gas], fs=1/self.stepSize, nperseg=windowsize, noverlap=windowsize-2, nfft=windowsize*5)
            t=t+self.windowSize*self.stepSize
            self.stftSBFre = f
            self.stftSBDelay = t
            vmax=abs(self.stftSB[gas][0:int(len(f)/10)]).max()*ratio
            vmin=abs(self.stftSB[gas]).min()*ratio
            #vmax=2*np.pi+0.5
            #vmin=-0.5
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            levels = np.arange(vmin,vmax,(vmax-vmin)/10)
            im=ax.contourf(t*10**15, f/10**12*33.35641,
                           np.abs(self.stftSB[gas]), levels=levels, cmap='jet', norm=norm)
            #im=ax.contourf(t*10**15, f/10**12*33.35641,
            #               np.angle(self.stftSB[gas])-np.angle(self.stftSB['rebin_Ch8']), levels=levels, cmap='jet')
            ax_histx.plot(t*10**15, np.sum(np.abs(self.stftSB[gas]),axis=0))
            ax_histx.set_xticklabels([])
            ax_histy.plot(np.mean(np.abs(self.stftSB[gas]),axis=1),f/10**12*33.35641)
            #ax_histy.plot(np.mean(np.angle(self.stftSB[gas])-np.angle(self.stftSB['rebin_Ch8']),axis=1),f/10**12*33.35641)
            #plt.xlim([100, self.delayB[-1]])
            ax_histy.set_ylim([200, 4500])
            ax.set_ylabel('Frequency [cm-1]')
            ax.set_xlabel('Time [fs]')
            #ax_histy.set_yticklabels([])
            ax_histy.yaxis.tick_right()
            #plt.clim(0,abs(self.stftSB[gas][:int(len(f)/15)]).max())
            #plt.colorbar(im, cax=a1, ticks=levels[::100])
            d = '_55_90'
            #plt.title(gas+d)
            fig.tight_layout()
            #plt.savefig("sfts_180fs.png",dpi=720, bbox_inches='tight', pad_inches=0.2)
            plt.show()

    def STFTS2(self, windowsize=300):
        startDelay = 90
        stepNum=0
        _preP = 'window'+'_'
        stft2 = {}
        while (startDelay+windowsize+(stepNum+1)*windowsize/2)*1E-15<self.delayB[-1]:
            stft2[str(startDelay+(stepNum+1)*windowsize/2)] = {}
            self.FFT3(windowSize=0, delayRange=[(startDelay+stepNum*windowsize/2)*1E-15,(startDelay+windowsize+stepNum*windowsize/2)*1E-15], rebinF=1,paddingF = 4.5, useWindow=True, zeroDirection='left', phaseCompensate=False, smooth=True,test = False)
            for gas in ['Ch8','Ch0','Ch2','Ch4','Ch6']:
                f_window = self.fftSB[_preP+'fre']
                Y_window = np.array([np.mean(np.abs(self.fftSB[_preP+gas+'_fft']),axis=0),np.std(np.abs(self.fftSB[_preP+gas+'_fft']),axis=0)])
                Y_window=np.where(f_window>=100,Y_window,0)/np.amax(Y_window)
                P_inter = np.angle(self.fftSB[_preP+gas+'_fft'])
                P_inter = np.mean(P_inter[:,np.where(np.logical_and(3656-30<f_window,f_window<3655+30))],2)
                P_inter = P_inter.flatten()
                P_inter=np.where(P_inter<0,P_inter+2*np.pi,P_inter)
                if gas == 'Ch8':
                    P_ref = P_inter
                    P_window =  np.array([np.mean(P_inter,axis=0), np.std(P_inter,axis=0)])
                else:
                    P_inter = P_inter-P_ref
                    P_window =  np.array([np.mean(P_inter,axis=0), np.sqrt(np.std(P_inter,axis=0)**2+np.std(P_ref,axis=0)**2)])
                P_window=P_window/np.pi
                #aa = len(f_window[(f_window<100)])
                #bb=len(f_window[(f_window<4150)])
                #f_window = f_window[aa:bb]
                #Y_window = Y_window[:,aa:bb]
                #P_window = P_window[:,aa:bb]
                #P_window[0]=np.where(P_window[0]<0,P_window[0]+2*np.pi,P_window[0])
                #P_window = P_window/np.pi
                stft2[str(startDelay+(stepNum+1)*windowsize/2)][gas]=P_window
                ##print(P_window)
                #if stepNum==0 and gas=='Ch8':
                #    stft2['fre'] = f_window
                #if gas == 'Ch2':
                #    plt.plot(f_window, self.baseLineRemove(Y_window[0]/np.amax(Y_window[0])))
                #    plt.errorbar(f_window,P_window[0],yerr=P_window[1], color='r', ecolor='r')
                #    plt.show()
            stepNum=stepNum+1
        color={}
        color['Ch0']='r'
        color['Ch2']='g'
        color['Ch4']='b'
        color['Ch6']='c'
        plt.figure(figsize=(7,5))
        for gas in ['Ch0','Ch2','Ch4','Ch6']:
            t=[]
            p=[]
            std=[]
            for key in stft2.keys():
                ##print(key)
                if key == 'fre':
                    continue
                t = t+[float(key)]
                if stft2[key][gas][0]+0.1<0:
                    p = p+[stft2[key][gas][0]+2]
                else:
                    p = p+[stft2[key][gas][0]]
                std=std+[stft2[key][gas][1]]
            plt.errorbar(t,p,yerr=std,ecolor='black',fmt=color[gas],mec=color[gas],mfc=color[gas],capsize=5, marker='s',markersize=10,elinewidth=2,markeredgewidth=2,label=self.label[gas])
            plt.xlabel('Delay(fs)')
            plt.ylabel('Phase($\pi$)')
            plt.legend()
        plt.tight_layout()
        plt.show()


    #def plotPhase(self):
    #    if ifsavePhase:
    #        self.wks = op.new_sheet('w',lname=str('Phase')+str('_')+self.folder)
    #    self.phase = {}
    #    ii=1
    #    for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8']:#
    #        self.phase[gas] = {}
    #        #if gas == 'Ch8':
    #        #    omega = [520.50412,619.22042,804.6874,879.47248,1136.73313,1334.16573,1588.43498,1687.15128,2144.83594,2318.33731,2668.33146,3200.80119,3643.52883]
    #        #elif gas == 'Ch0':
    #        #    omega = [520.50412,804.6874,1334.16573,1588.43498,2144.83594,3643.52883]
    #        #elif gas == 'Ch2' or gas == 'Ch4':
    #        #    omega = [3643.52883]
    #        #elif gas == 'Ch6':
    #        #    omega = [520.50412,804.6874,1334.16573,1588.43498,2144.83594,3643.52883]
    #        omega = [526.49165,625.20883,699.99458,810.67748,882.47179,1145.71762,1343.15199,1594.43209,1696.1407,2153.82946,2324.34096,2677.32968,3212.79562,3655.52723]
    #        mean = []
    #        std = []
    #        for f in omega:#526.49165,625.20883,699.99458,810.67748,882.47179,1145.71762,1343.15199,1594.43209,1696.1407,2153.82946,2324.34096,2677.32968,3212.79562,3655.52723
    #            _inter = []
    #            for i in self.phaseBottle[gas].keys():
    #                #if f==3643.52883:
    #                #    plt.plot(self.phaseBottle[gas][str(i)][str(f)])
    #                #    plt.show()
    #                #A, w, p, c = self.phaseBottle[gas][str(i)][str(f)+'_fit'][1]
    #                #if f==3643.52883:
    #                #    plt.plot(self.phaseBottle[gas][str(i)][str(f)][:int(33356.40952/f/self.delayStep/1e15)],label=str(i))
    #                #p = sps.argrelextrema(self.phaseBottle[gas][str(i)][str(f)], np.greater)[0][0]/(33356.40952/f/self.delayStep/1e15)*np.pi*2
    #                p = np.argmax(self.phaseBottle[gas][str(i)][str(f)][:int(33356.40952/f/self.stepSize/1e15)])/int(33356.40952/f/self.stepSize/1e15)*np.pi*2
    #                if i =='9' or i == '10':
    #                    continue
    #                if np.mean(_inter)-p>5:
    #                    _inter = _inter + [p+np.pi*2]
    #                elif np.mean(_inter)-p<-5:
    #                    _inter = _inter + [p-np.pi*2]
    #                else:
    #                    _inter = _inter + [p]
    #            #plt.plot(_inter)
    #            #plt.show()
    #            _std = np.std(_inter)
    #            if _std>1:
    #                #mean = mean + ['']
    #                #std = std + ['']
    #                _mean = np.mean(_inter)
    #                mean = mean + [_mean]
    #                std = std + [_std]
    #            else:
    #                _mean = np.mean(_inter)
    #                mean = mean + [_mean]
    #                std = std + [_std]
    #        self.phase[gas] = [omega, mean, std]
    #        if ifsavePhase:
    #            self.wks.from_list(0, omega, 'X')
    #            self.wks.from_list(ii, mean, lname=self.label[gas], axis='Y')
    #            self.wks.from_list(ii+1, std, lname=self.label[gas], axis='E')
    #            ii=ii+2
    #        plt.errorbar(omega,np.array(mean)%(2*np.pi),yerr=np.array(std), fmt='s',label = self.label[gas])
    #    plt.legend()
    #    plt.ylim([0,6.28])
    #    #plt.savefig(os.path.join(os.path.join(self.savePath,r'phase.png')),dpi=720,bbox_inches='tight',pad_inches=0,transparent=True)
    #    plt.show()

    def cal_ratio(self):
        self.ratio = {}
        #if volts out=1V, volts in: 
        #1V, 10mV, 10mV, 1V ,5V and 10mV #4.5E+14_H2O 9e-7 mbar
        #1V, 20mV, 10mV, 1V ,5V and 20mV #7.2E+14_H2O 5e-7 mbar
        #1V, 20mV, 10mV, 2V ,5V and 20mV #8.9E+14_H2O 5e-7 mbar
        boxcarA=[2,5]
        u = np.sum(self.specBigBottleB['Ch0'][-9000:])*2

        #print(u)
        i=0
        for gas in ['Ch0', 'Ch2']:
            self.ratio[gas] = np.sum(self.specBigBottleB[gas][-9000:])/u*boxcarA[i]
            i=i+1

    def fit_phase(self, omega):
        '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
        dw=self.dw
        def fitFunc(params, x, data, omega):
            amp = params['amp'].value
            phi = params['phi'].value
            w0 = params['w0'].value
            n = params['n'].value
            
            model = amp*np.sinc(np.pi*(x-w0)/dw/n)*np.exp(1j*(np.pi*(x-w0)/dw/n+phi))
            #model = amp*np.sin(np.pi*self.L*(x-w0))/(2*np.pi*self.L*(x-w0)*(1-self.L*self.L*(x-w0)*(x-w0)))*np.exp(1j*(np.pi*(x-w0)/self.dw/n+phi))
            ##print('here self.dw is ' +str(self.dw))
            extre = sps.argrelextrema(x, np.greater)
            weight = np.zeros(np.size(x))+0.1
            for i in extre:
                if np.abs(x[i]-omega)<50:
                    weight[i-1]=weight[i-1]+100
                    weight[i]=weight[i]+100
                    weight[i+1]=weight[i+1]+100
                else:
                    weight[i-1]=weight[i-1]+10
                    weight[i]=weight[i]+10
                    weight[i+1]=weight[i+1]+10
            res=np.abs(np.real((model-data))*weight)+1j*np.abs(np.imag((model-data)*weight))
            return res.view(float)#that's what you want to minimize
        # create a set of Parameters
        if omega == 3657:
            params = Parameters()
            params.add('amp', value= 45.5, min=1,max=100 ,brute_step=0.00000001) #value is the initial condition
            params.add('phi', value=1, min=0, max=3*np.pi ,brute_step=0.00000001)#value=-2
            params.add('w0', value= 3655, min=3650, max=3660 ,brute_step=0.0000001)
            params.add('n', value= 1.0, min=1.07, max=1.0701 ,brute_step=0.0000001)
            

        elif omega == 3221:
            params = Parameters()
            params.add('amp', value= 2e5, min=1,max=1e6 ,brute_step=0.001) #value is the initial condition
            params.add('phi', value=1, min=0, max=3*np.pi ,brute_step=0.001)
            params.add('w0', value= 3216, min=3100, max=3300 ,brute_step=0.001)#3217,3214,3217
            params.add('n', value= 1.032, min=1.032, max=1.033 ,brute_step=0.001)
            
        elif omega ==2152:
            params = Parameters()
            params.add('amp', value= 2e4, min=1,max=1e6 ,brute_step=0.001) #value is the initial condition
            params.add('phi', value=1.5, min=0, max=3*np.pi, brute_step=0.001)
            params.add('w0', value= 2151, min=2100, max=2200 ,brute_step=0.001)
            params.add('n', value= 1, min=0.9, max=1.2 ,brute_step=0.001)

        elif omega == 1342:
            params = Parameters()
            params.add('amp', value= 10, min=1,max=100 ,brute_step=0.001) #value is the initial condition
            params.add('phi', value=0.6, min=0, max=3*np.pi ,brute_step=0.001)
            params.add('w0', value= 1344, min=1300, max=1400 ,brute_step=0.001)
            params.add('n', value= 1, min=0.9, max=1.2 ,brute_step=0.001)
            
        elif omega ==810:
            params = Parameters()
            params.add('amp', value= 2e4, min=1,max=1e6 ,brute_step=0.001) #value is the initial condition
            params.add('phi', value=1.5, min=0, max=2*np.pi, brute_step=0.001)
            params.add('w0', value= 811.8, min=750, max=850 ,brute_step=0.001)
            params.add('n', value= 1, min=0.9, max=1.2 ,brute_step=0.001)

        elif omega == 4161:
            params = Parameters()
            params.add('amp', value= 45.5, min=1,max=100 ,brute_step=0.00000001) #value is the initial condition
            params.add('phi', value=4.8, min=0, max=2*np.pi ,brute_step=0.00000001)#value=-2
            params.add('w0', value= 4161.5, min=4100, max=4200 ,brute_step=0.0000001)
            params.add('n', value= 1.0, min=1.07, max=1.0701 ,brute_step=0.0000001)           

        #######################################
        self.interfftSB = load_obj(os.path.abspath(os.path.join(d.savePath, str(self.smooth)+str(self.folder)+r'_fftSB.pkl')))
        _preP = 'window'+'_'
        gas = 'Ch8'
        f = self.fftSB[_preP+'fre']
        Y = self.fftSB[_preP+gas+'_fft']

        
        aa = len(f[(f<omega-150)])
        bb=len(f[(f<omega+150)])
        #aa = len(f_window[(f_window<100)])
        #bb=len(f_window[(f_window<50000)])
        f_window = f[aa:bb]
        Y_window = Y[:,aa:bb]
        Y_window_re = np.real(Y_window)
        Y_window_im = np.imag(Y_window)
        delta = 1
        x=f_window#[delta:]
        ###########################################
        phase = []
        for i in range(np.shape(Y)[0]):
        
            #y=Y_window_re[i,:-delta]+1j*Y_window_im[i,delta:]
            y=Y_window_re[i]+1j*Y_window_im[i]
            y = y/np.max(np.abs(y))
            # do fit, here with leastsq model
            result = minimize(fitFunc, params, args=(x,y,omega),method='least_squares')

            # calculate final result result.params['amp'].value
            final = result.params['amp'].value*np.sinc(np.pi*(x-result.params['w0'].value)/dw/result.params['n'].value
                                                       )*np.exp(1j*(np.pi*(x-result.params['w0'].value)/dw/result.params['n'].value+result.params['phi'].value))
            # write error report
            result.params.pretty_print()
            phase = phase + [result.params['phi'].value]
            
            plt.plot(x, np.real(y),'r', label = 'Data')
            plt.plot(x, np.real(final), 'b',label = 'Fitting')
            plt.plot(x, np.imag(y),'r--',label = 'Data')
            plt.plot(x, np.imag(final), 'b--',label = 'Fitting')
            plt.plot(x, np.abs(y), 'g')
            plt.legend()
            plt.show()
        #print(phase)
        phase = np.array(phase)/np.pi#+self.windowSize/self.inputLength*omega/self.dw*2/dw
        #print(self.windowSize)
        #print(self.windowSize/self.inputLength*omega/self.dw*2)
        print([omega, np.mean(phase)%2, np.std(phase)])
        #np.savetxt(str(omega)+'.txt',np.stack((x,y),axis=1))
        np.savetxt(str(omega)+'_fitParameter.txt',[omega, np.mean(phase)%2, np.std(phase)])
        return phase

if __name__ == '__main__':
    for ff in [r'3.6E+14_H2',r'6.0E+14_H2',r'8.0E+14_H2']:#r'3.6E+14_H2',r'6.0E+14_H2',r'8.0E+14_H2'True
        d = FFT_ionS(ff)
        if d.checkSavedData():
            d.read()
            d.delayCorrection()
        d.transition()
        d.cal_ratio()
        print(d.ratio)
        d.findZeroDelay3()
        d.show_Spectra(ifsaveT=False)
        #d.FFT3(windowSize=100, delayRange=[300*1E-15,1000*1E-15], rebinF=1,paddingF = 5, useWindow=True, zeroDirection='left', phaseCompensate=False, smooth=True,test = False)
        #d.FFT3(windowSize=90, delayRange=False, rebinF=1,paddingF = 20, useWindow=True, zeroDirection='left', phaseCompensate=True, smooth=True,test = False)
        #d.show_FFT(ifsaveFFT=False,ifsaveT=False)
        #