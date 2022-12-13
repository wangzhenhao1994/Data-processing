from fileinput import filename
from weakref import ref
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from yaml import load
#mpl.style.use('seaborn-poster')
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
import pickle
from decimal import Decimal
from cal_intensity import cal_intensity
from calculate_k_b import Calibration_mass
from BaselineRemoval import BaselineRemoval
from lmfit import minimize, Parameters, Parameter, report_fit
import originpro as op

ifsave = False
ifsaveT = False
ifsavePhase = False
if ifsave or ifsaveT or ifsavePhase:
    op.set_show(show=True)
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
        self.rootPath= pl.PureWindowsPath(r'D:\DataProcessing\CH2O')
        #self.savePath= pl.PureWindowsPath(r'C:\Users\user\Desktop\Data_newTOF\dataProcessing\4.5E+14_H2O')
        #self.savePath= pl.PureWindowsPath(os.path.join(r'C:\Users\user\Desktop\Data_newTOF\dataProcessing\09122022\H2O',folder))
        self.savePath= pl.PureWindowsPath(os.path.join(r'D:\DataProcessing\CH2O',folder))
        self.delayB,self.stepSize = np.linspace(start=0,stop=1300, num=13000,endpoint=False,retstep=True)
        self.delayB = self.delayB*10**-15
        self.stepSize = self.stepSize*10**-15
        self.delayB_noPad = np.linspace(start=0,stop=1300, num=13000)*10**-15
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
        self.molecule = str('H20')

        #setting for data from digitizer
        self.channelSize = 12032#24000#1536
        self.scanLengthD = 3320#1200
        self.peakRange = [-100, 100]  # range of the peak
        self.delayD = np.arange(self.scanLengthD)/self.scanLengthD*100*2*2*3.33564*10**-15#*3602/3647
        self.spectraBottleD = {}
        self.fftSD = {}
        self.stftSD = {}
        self.dataD = 0

    def read(self):
        totalCounter = {}
        #for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
        #    self.fftSB[gas] = 0
        for gas in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
            self.specBigBottleB[gas] = np.zeros((1000,13000))
            totalCounter[gas] = 0
        for file in os.listdir(self.savePath):
            filename = os.fsdecode(file)
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

    def window(self, windowSize=0, direction='left', useWindow = True):
        '''
        windowSize is in fs.
        '''
        #windowSize = windowSize+np.abs(self.delayB[0])*1e15 #????????????????????????
        for gas in list(self.specBigBottleB):
            self.specBigBottleB[gas], self.specBigBottleB['window_'+gas], self.delayB = self.inter_window(self.specBigBottleB[gas], self.delayB, windowSize=windowSize, direction=direction, useWindow=useWindow)
           
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
        data = data[:int(np.size(data)/2)*2]
        __len = np.size(data)
        windowSize = int(windowSize*1e-15/self.stepSize)
        self.windowSize = windowSize
        if smooth:
            data=sps.savgol_filter(data, window_length=10, polyorder=3,deriv=2,delta=2, mode='interp')
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
        window=polynomial(window, order=1, plot=False)
        window = window - window[-1]#shift the base line to zero

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
        for i in range(self.specBigBottleB['Ch8'].shape[0]):
            _Spec = self.specBigBottleB['Ch8'][i]
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

    def calDrift(self, _cRange):
        plt.clf()
        _iS = np.array(np.zeros(self.specBigBottleB['Ch0'].shape))
        sumSpec=_iS
        for gas in self.specBigBottleB.keys():
            sumSpec=sumSpec+self.specBigBottleB[gas]

        for i in range(self.specBigBottleB['Ch0'].shape[0]):
            #_iS[i]=self.butter_bandpass_filter(self.specBigBottleB[gas][i], (3643.52883-0.1)/33.35641*1e12, (3643.52883+0.1)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            _iS[i]=self.butter_bandpass_filter(sumSpec[i], 500/33.35641*1e12, 2000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))

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
            print('Data is too little!')
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
            if _shift.size <4:
                empty = np.zeros((3,self.specBigBottleB[gas].shape[1]))
                self.specBigBottleB[gas] = np.concatenate((self.specBigBottleB[gas], empty), axis=0)
            iS = sp.interpolate.RectBivariateSpline(range(self.specBigBottleB[gas].shape[0]), self.delayB_noPad*1e15, self.specBigBottleB[gas])
            inter = np.zeros((oShape,self.specBigBottleB[gas].shape[1]))
            for i in range(oShape):
                inter[i] = iS.ev(i, self.delayB_noPad*1e15+_shift[i])
            iinter[gas] = inter
        self.specBigBottleB = iinter

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
            #vmax=abs(self.stftSB[gas][0:int(len(f)/10)]).max()*ratio
            #vmin=abs(self.stftSB[gas]).min()*ratio
            vmax=2*np.pi+0.5
            vmin=-0.5
            #norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            levels = np.arange(vmin,vmax,(vmax-vmin)/10)
            #im=ax.contourf(t*10**15, f/10**12*33.35641,
            #               np.abs(self.stftSB[gas]), levels=levels, cmap='jet', norm=norm)
            im=ax.contourf(t*10**15, f/10**12*33.35641,
                           np.angle(self.stftSB[gas])-np.angle(self.stftSB['rebin_Ch8']), levels=levels, cmap='jet')
            ax_histx.plot(t*10**15, np.sum(np.abs(self.stftSB[gas]),axis=0))
            ax_histx.set_xticklabels([])
            #ax_histy.plot(np.mean(np.abs(self.stftSB[gas]),axis=1),f/10**12*33.35641)
            ax_histy.plot(np.mean(np.angle(self.stftSB[gas])-np.angle(self.stftSB['rebin_Ch8']),axis=1),f/10**12*33.35641)
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
                #print(P_window)
                #if stepNum==0 and gas=='Ch8':
                #    stft2['fre'] = f_window
                #if gas == 'Ch2':
                #    plt.plot(f_window, self.baseLineRemove(Y_window[0]/np.amax(Y_window[0])))
                #    plt.errorbar(f_window,P_window[0],yerr=P_window[1], color='r', ecolor='r')
                #    plt.show()
            stepNum=stepNum+1

        for gas in ['Ch0','Ch2','Ch4','Ch6']:
            t=[]
            p=[]
            std=[]
            for key in stft2.keys():
                #print(key)
                if key == 'fre':
                    continue
                t = t+[float(key)]
                if stft2[key][gas][0]+0.1<0:
                    p = p+[stft2[key][gas][0]+2]
                else:
                    p = p+[stft2[key][gas][0]]
                std=std+[stft2[key][gas][1]]
            plt.errorbar(t,p,yerr=std,ecolor='black',fmt='r',capsize=5, marker='s',markersize=10,mec='r',mfc='r',elinewidth=2,markeredgewidth=2,)
            plt.xlabel('Delay(fs)')
            plt.ylabel('Phase($\pi$)')
            plt.show()
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
            print(gas)
            n=1
            self.phaseSpecBottleB[gas] = self.specBigBottleB[gas][-int(self.specBigBottleB[gas].shape[0]/n)*n:].reshape(int(self.specBigBottleB[gas].shape[0]/n), n, self.specBigBottleB[gas].shape[1]).sum(axis=1)
            self.specBigBottleB[gas] = self.specBigBottleB[gas].sum(axis = 0)#np.take(self.specBigBottleB[gas],[i for i in range(self.specBigBottleB[gas].shape[0]) if i not in self.goodSpecIndex],axis=0).sum(axis = 0)
        #self.specBigBottleB['Ch4'] = self.specBigBottleB['Ch4']-self.specBigBottleB['Ch10']
        #self.specBigBottleB['Ch6'] = self.specBigBottleB['Ch6']-self.specBigBottleB['Ch0']
        
        self.label = {}
        self.label['Ch0'] = 'Mass1'
        self.label['Ch2'] = 'Mass2'
        self.label['Ch4'] = 'Mass28'
        self.label['Ch6'] = 'Mass29'
        self.label['Ch8'] = 'Mass30'  
 
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
        self.stepSize = self.stepSize*rebinF
        _size = self.specBigBottleB['Ch0'].size
        paddingSize = int(_size*paddingF)
        self.paddingSize = paddingSize
        for gas in self.phaseSpecBottleB.keys():
            for i in range(self.phaseSpecBottleB[gas].shape[0]):
                interSpec = self.phaseSpecBottleB[gas][i][-self.specBigBottleB[gas].size:]
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
    
    
    def show_FFT(self):
        self.dcRange=0
        fig, ax = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True,figsize=(18,9))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel("Frequency ($\mathrm{cm^{-1}}$)")
        plt.ylabel("Amplitude (a.u.)", labelpad=25)
        plt.text(1.08,0.38,'Relative Phase v.s. Mass 18 ($\pi$)', rotation = 270)
        plt.title(self.folder,fontsize=40,pad =20)
        axPP = []
        for axx in ax.flatten():
            axPP = axPP + [axx.twinx()]
        axPP[0].get_shared_y_axes().join(axPP[0], axPP[1],axPP[2])

        lab=['Mass 18','Mass 1','Mass 2','Mass 16','Mass 17']
        i=0
        if ifsave:
            self.wks = op.new_sheet('w',lname=str('FFT')+str('_')+self.folder)

        _preP = 'window'+'_'
        for gas in ['Ch8','Ch0','Ch2','Ch4','Ch6']:
            if 'filter' in gas  or 'window' in gas or 'rebin' in gas or 'com' in gas:
                continue
            elif gas not in ['Ch8','Ch0','Ch2','Ch4','Ch6']:
                continue
            
            axF=ax[int(i/3),i%3]
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
            if gas == 'Ch8':
                P_ref = P_inter
                P_window =  np.array([np.mean(P_inter,axis=0), np.std(P_inter,axis=0)])
            else:
                P_inter = P_inter-P_ref
                P_window =  np.array([np.mean(P_inter,axis=0), np.sqrt(np.std(P_inter,axis=0)**2+np.std(P_ref,axis=0)**2)])

            
            #aa = len(f_window[(f_window<100)])
            #bb=len(f_window[(f_window<4000)])
            aa = len(f_window[(f_window<200)])
            bb=len(f_window[(f_window<4500)])
            f_window = f_window[aa:bb]
            Y = Y[:,aa:bb]
            Y_window = Y_window[:,aa:bb]
            Y_window[0] = self.baseLineRemove(Y_window[0]/np.amax(Y_window[0]))
            Y_window_im = Y_window_im[:,aa:bb]
            Y_window_re = Y_window_re[:,aa:bb]
            P_window = P_window[:,aa:bb]
            P_window[0]=np.where(P_window[0]<0,P_window[0]+2*np.pi,P_window[0])
            P_window = P_window/np.pi
            P_window[0] = np.where(Y_window[0]>0.05,P_window[0],np.inf)
            #fitRes = self.fit_phase(f_window,Y[0]/np.amax(np.abs(Y[0])),[3655.52723])
            axF.plot(f_window, Y_window[0], label=self.label[gas])
            #axF.plot(f_window, Y_window_re[0]/(np.amax(Y_window_re[0])-np.amin(Y_window_re[0])), label=label+'_re')
            #axF.plot(f_window, Y_window_im[0]/(np.amax(Y_window_im[0])-np.amin(Y_window_im[0])), label=label+'_im')
            #axF.plot(f_window, self.baseLineRemove(Y_window_re[0]/np.amax(Y_window[0])))#, label=label+'_re')
            #axF.plot(f_window, self.baseLineRemove(Y_window_im[0]/np.amax(Y_window[0])))#, label=label+'_im')
            if gas == 'Ch8':
                peaks, _ = sps.find_peaks(Y_window[0],width=5,prominence=0,distance=40,height=0.05)
                axF.plot(f_window[peaks],Y_window[0][peaks],'+r')
                for peak in peaks:
                    axF.text(f_window[peak],Y_window[0][peak]+0.05,str(int(f_window[peak])))
                pass
            else:
                axP.errorbar(f_window,P_window[0],yerr=P_window[1], color='r', ecolor='r')
            #plot(f_window,P_window,'r')
            #axF.set_ylim([0,1])
            #axF.set_xlim([200,4500])
            axP.set_ylim([-0.5,2.5])
            axF.legend(loc=1,ncol=2,shadow=False,framealpha=0,scatterpoints=0)
            i=i+1
            if ifsave:
                self.wks.from_list(0, f_window, lname="Frequency", axis='X')
                self.wks.from_list(i, np.abs(Y_window), lname=label, axis='Y')

        fig.tight_layout()
        #plt.savefig(os.path.join(os.path.join(self.savePath,self.folder+'.png')),dpi=720,bbox_inches='tight',pad_inches=0,transparent=True)
        plt.show()

    def show_Spectra(self, shift=0):
        gs = gridspec.GridSpec(2, 3)
        #gs = gridspec.GridSpec(1, 1)
        fig = plt.figure(figsize=(20,8))
        #ax = fig.add_subplot(111)
        lab=['Mass1','Mass2','Mass16','Mass17','Mass18','Ch2+Ch4','Mass1+Mass16', 'Mass1+Mass17', 'Mass16+Mass17','Mass17+Mass18']
        i=0
        if ifsaveT:
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
            ax.plot(delay*10**15+shift,
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
                self.wks.from_list(i, self.specBigBottleB['filter_'+gas], lname=gas, axis='Y')
        #plt.legend()
        fig.tight_layout()
        plt.show()

if __name__ == '__main__':
    directory = os.path.abspath(r'D:\DataProcessing\CH2O')
    for ff in [name for name in os.listdir(directory)
            if os.path.isdir(os.path.join(directory, name))]:
        print(ff)
    #for ff in [r'pu1.1E+15pr1.1E+15_CH2O']:
    #for ff in [r'pu9.2E+13pr5.3E+14_CH2O',r'pu9.2E+13pr9.1E+14_CH2O',r'pu9.2E+13pr1.4E+15_CH2O']:#benging?
    #for ff in [r'pu1.2E+15pr2.2E+14_CH2O',r'pu1.2E+15pr8.7E+14_CH2O']:
        d = FFT_ionS(ff)
        if d.checkSavedData():
            d.read()
            d.delayCorrection()
        d.transition()
        d.findZeroDelay3()
        #d.show_Spectra()
        d.FFT3(windowSize=90, delayRange=False, rebinF=1,paddingF = 5, useWindow=True, zeroDirection='left', phaseCompensate=False, smooth=False,test = False)
        d.show_FFT()
        #mdic = {"Ch8": d.specBigBottleB['Ch8'], "Ch0": d.specBigBottleB['Ch0'],"label": "experiment"}
        #from scipy.io import savemat
        #savemat("matlab_matrix.mat", mdic)
        #d.rebinS(5)
        #d.window(windowSize=200,useWindow=False)
        #d.STFTS(gas='rebin_Ch8',windowsize=400)
        #d.STFTS(gas='rebin_Ch0',windowsize=400)
        #d.STFTS2(windowsize=182.5*2)
