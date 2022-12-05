from fileinput import filename
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
from Superlets.python.superlet import *

import originpro as op
ifsave = False
ifsaveT = False
ifsavePhase = False
if ifsave or ifsaveT or ifsavePhase:
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
    #def __init__(self, filename, massRange=[5, 90], dcRange=2, cal_mass=[45,61], cal_pixel=[415,1141]):
    def __init__(self, filename, scanTime, sampRate, massRange=[5, 50], molecule='H2O', intensity = 0, dcRange=5, ChN=11, cal_mass=[2,16], cal_pixel=[3082,10216]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[494,583]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[197,279]):
        '''

        '''
        self.filename = filename
        self.saveRef = str(filename)
        self.filepath = []
        self.rootPath= pl.PureWindowsPath(r'C:\Users\user\Desktop\Data_newTOF')
        self.scanTime = scanTime
        self.trueScanTime = ceil(self.scanTime/1.6)*1.6
        self.scanLengthB = 3360
        self.delayB = np.arange(self.scanLengthB)/self.scanLengthB*100*2*2*3.33564*10**-15*3642/3647
        self.rebin_delay = None
        self.ifrebin = False
        self.longstage = True
        self.data = 0
        self.massRange = massRange
        self.dcRange = dcRange
        self.spectraBottleB = {}
        self.spectraBottleB_noPad = {}
        self.interSpectraBottleB = {}
        self.correctedInterSpecB = {}
        self.fftSB = {}
        self.stftSB = {}
        self.waveletsCoe = {}
        self.zeroIndex = {}
        self.zeroIndexMin = 10000000000000000000000000
        self.zeroIndexMax = 0
        self.retrivedPhase = {}
        self.spec = None
        self.scanNo = None
        self.ChN = ChN
        self.molecule = molecule
        try:
            self.intensity = 'pu'+str('%.1E' % Decimal(intensity[0]))+'pr'+str('%.1E' % Decimal(intensity[1]))
        except TypeError:
            self.intensity = str('%.1E' % Decimal(intensity))  
        self.stage = 'piezo'
        self.savePath = os.path.join(self.rootPath, r'dataProcessing', str(self.intensity)+'_'+self.molecule)

        #setting for data from digitizer
        self.calculator = Calibration_mass(cal_mass, cal_pixel)
        self.channelSize = 12032#24000#1536
        self.scanLengthD = 3320#1200
        self.peakRange = [-100, 100]  # range of the peak
        self.delayD = np.arange(self.scanLengthD)/self.scanLengthD*100*2*2*3.33564*10**-15*3602/3647
        self.gasBottle = {
            "O+": self.calculator.cal_pixel(16)+self.peakRange,
            "H2+": self.calculator.cal_pixel(2)+self.peakRange,
            "H+": self.calculator.cal_pixel(0.99)+self.peakRange,
            "O+H2": #[0,self.channelSize]
            np.append(self.calculator.cal_pixel(16) + self.peakRange, self.calculator.cal_pixel(2) + self.peakRange),
            "O-H2": #[0,self.channelSize]
            np.append(np.append(self.calculator.cal_pixel(16) + self.peakRange, self.calculator.cal_pixel(2) + self.peakRange), 0),
        }
        self.spectraBottleD = {}
        self.fftSD = {}
        self.stftSD = {}
        self.dataD = 0

        
    def pathFinder(self):
        for fn in self.filename:
            interpath = pl.PureWindowsPath(fn[9:16], fn[9:19])
            self.filepath = self.filepath + [pl.PureWindowsPath(self.rootPath, interpath, fn + r'.hdf5')]
            #self.exportPath = pl.PureWindowsPath(self.rootPath, r'plot20220602')
            #self.exportPath = pl.PureWindowsPath(self.rootPath, interpath, r'processed')
        #if not os.path.exists(self.exportPath):
        #    os.mkdir(self.exportPath)

    def checkData(self):
        #with h5py.File(self.filename, 'r+') as f:
        #    print(np.array(f['parameters']))
        for fp in self.filepath:
            try:
                with h5py.File(fp, 'r+') as f:
                    bigger = True
                    i = 0
                    while bigger:
                        try:
                            a = f['dataD'][i*self.channelSize,0]
                            i+=1
                        except ValueError:
                            bigger = False
                            self.scanNo = i
                            try:
                                f.create_dataset('scanNo', data=i)
                            except RuntimeError:
                                del f['scanNo']
                                f.create_dataset('scanNo', data=i)
                            print('The number of scan is '+str(i)+'!')
            except OSError:
                pass
    
    def combineData(self):
        pass

    def read_splitB(self, isCorrection = True, sumN=1):
        ChN=self.ChN #number of channels measured
        for fp in self.filepath:
            self.delayB = np.arange(self.scanLengthB)/self.scanLengthB*100*2*2*3.33564*10**-15*3642/3647
            if isCorrection:
                with h5py.File(fp, 'r+') as f:
                    print(f.keys())
                    data = np.array(f['dataB'], dtype=float)
                    print(data.shape)
                    m = int(data.shape[0]/self.ChN/sumN)
                    for i in [0,2,4,6,8,10]:
                        self.interSpectraBottleB['Ch'+str(i)] = data[i::ChN][-int(m*sumN):].reshape(m, sumN, data.shape[1]).sum(axis=1)
                    self.interinterSpectraBottleB = {}
                    num = int(m/59)
                    for gas in self.interSpectraBottleB.keys():
                        self.interinterSpectraBottleB[gas] = {}
                        self.spectraBottleB[gas] = np.zeros((num, 13000))
                    for i in range(num):
                        for gas in self.interSpectraBottleB.keys():
                            _inter = self.interSpectraBottleB[gas][i*59:(i+1)*59]
                            self.interinterSpectraBottleB[gas][str(i)] = _inter   
            self.num = num
            self.delayCorrection()
            if not os.path.exists(self.savePath):
                os.mkdir(self.savePath)
            print(self.interSpectraBottleB['Ch0'].shape)
            save_obj(self.interSpectraBottleB, pl.PureWindowsPath(self.savePath, os.path.split(fp)[1].replace(r'.hdf5','')+r'.pkl'))

    def read_splitD(self, overNight = True, firstTry = False, sumNo = 89, usefulRange = [0,2], cu = False):
        for fp in self.filepath:
            with h5py.File(fp, 'r+') as f:
                print(np.array(f['dataD'].shape))
                print(f.keys())
                dataSet_created={}
                if ('scanNo' in f):
                    self.scanNo = np.array(f['scanNo'], dtype=int)
                else:
                    self.checkData()
                if overNight:
                    splitNo = int(self.scanNo/sumNo)
                    for i in range(splitNo):
                        dataSet_created['sum'+str(i)]=[i*sumNo, (i+1)*sumNo]
                else:
                    dataSet_created = {'sum0':[0,self.scanNo-1]}
                for s in dataSet_created.keys():
                    e = s in f
                    if e and (not firstTry):
                        print('dataset \''+ s +'\' exists!!!')
                        pass
                    else:

                        print('Begin to create dataset \''+ s +'\'!')
                        data = np.zeros((self.channelSize, self.scanLengthD))
                        for i in range(dataSet_created[s][1],dataSet_created[s][0],-1):
                            try:
                                data = data + \
                                    np.array(f['dataD'][i*self.channelSize:(i+1)*self.channelSize])
                            except ValueError:
                                break
                        try:
                            f.create_dataset(s, data=data, compression="gzip") #create dataset 'sum'
                            print('Dataset \''+ s +'\' is created!')
                        except RuntimeError:
                            print('Dataset \''+ s +'\''+' already exists!')
                            try:
                                del f[s]
                                f[s] = data
                                print('Dataset \'' + s + '\''+' is overwrited!')
                            except RuntimeError:
                                print('Fail to overwrite dataset \'' + s + '\' !')
                if overNight and 'useful' in f and not firstTry and not cu:
                    self.data = np.array(f['useful'])
                    print(np.array(f['useful_description']))
                else:
                    dataSet_use = ['sum'+str(i) for i in range(usefulRange[0],usefulRange[1])]#[x for x in dataSet_created.keys()]
                    print('Now the dataset \'',  dataSet_use, '\' is in use!')
                    for d in dataSet_use:
                        self.data = self.data + np.array(f[d])
                    try:
                        f.create_dataset('useful', data=self.data) #create dataset 'sum'
                        useful_description = 'Sum of trace in range '+str((np.array(usefulRange)-[0,1])*sumNo)+' is in use!'
                        print(useful_description)
                        f.create_dataset('useful_description', data=useful_description) #create dataset 'sum'
                        print('Dataset \''+ 'useful' +'\' is created!')
                    except RuntimeError:
                        print('Dataset \''+ 'useful' +'\''+' already exists!')
                        if cu:
                            del f['useful']
                            print('Begin to overwrite dataset \''+ 'useful' +'\''+'!')
                            f['useful'] = self.data
                #if np.sum(self.dataD) == 0:
                #self.dataD = self.data
                #else:
                #    pass


    def mass_spectra(self):
        zeroIndex = self.findZeroDelayD()
        zeroIndexS=zeroIndex
        zeroIndexS.sort()
        zeroIndexMin = zeroIndexS[0]
        zeroIndexMax = zeroIndexS[-1]
        fn=0
        for fp in self.filepath:
            with h5py.File(fp, 'r+') as f:
                self.dataD = self.dataD + np.array(f['useful'])[:,(zeroIndex[fn]-zeroIndexMin):(self.scanLengthD-(zeroIndexMax-zeroIndex[fn]))]
            fn+=1
        for gas in self.gasBottle.keys():
            if len(self.gasBottle[gas])==2:
                [pixelMin, pixelMax] = list(map(int, self.gasBottle[gas]))
                #print(pixelMin, pixelMax)
                inter= np.sum(
                    self.dataD[pixelMin:pixelMax, :], 0)
                #print(pixelMin, pixelMax)
                self.spectraBottleD[gas] = inter
            elif len(self.gasBottle[gas])==4:
                [pixel0, pixel1, pixel2, pixel3] = list(map(int, self.gasBottle[gas]))
                inter1=np.sum(self.dataD[pixel0:pixel1, :], 0)
                inter2=np.sum(self.dataD[pixel2:pixel3, :], 0)
                self.spectraBottleD[gas] =inter1 + inter2
                print('Sum two mass!!')
            elif len(self.gasBottle[gas])==5:
                [pixel0, pixel1, pixel2, pixel3, useless] = list(map(int, self.gasBottle[gas]))
                inter1=np.sum(self.dataD[pixel0:pixel1, :], 0)
                inter2=np.sum(self.dataD[pixel2:pixel3, :], 0)
                self.spectraBottleD[gas] =inter1 - inter2
                print('Sub two mass!!')

    def FFT(self, t, y): #https://www.researchgate.net/post/What-are-the-basic-differences-between-FFT-and-DFT-and-DCT
        
        n = len(t)
        #y=np.roll(y,int(n/2)) #circular shift
        delta = (max(t) - min(t)) / (n)
        f = sft.fftshift(sft.fftfreq(n, delta))/ 10**12  # frequency unit THz
        fft_y = sft.fftshift(sft.fft(y))
        fft_y = fft_y/np.mean(np.abs(fft_y))
        Y = np.abs(fft_y)
        Y_r = np.real(fft_y)
        Y2 = np.where(Y<np.mean(Y)*0,0,fft_y)
        #plt.plot(np.abs(Y2))
        #plt.show()
        P = np.angle(Y2)
        return np.array([f, Y, P])

    def window(self, windowSize=0, direction='left'):
        '''
        windowSize is in fs.
        '''
        self.spectraBottleB_noPad=self.spectraBottleB.copy()
        self.delayB_noPad = self.delayB.copy()
        windowSize = windowSize+np.abs(self.delayB[0])*1e15
        self.windowSize = windowSize
        for gas in list(self.spectraBottleB):
            if 'rebin' not in gas:
                self.spectraBottleB[gas], self.spectraBottleB['window_'+gas], self.delayB = self.inter_window(self.spectraBottleB[gas], self.delayB, windowSize=windowSize, direction=direction)
            else:
                self.spectraBottleB[gas], self.spectraBottleB['window_rebin_'+gas], self.rebin_delay = self.inter_window(self.spectraBottleB[gas], self.rebin_delay, windowSize=windowSize, direction=direction)
        for gas in list(self.spectraBottleD):
            if 'rebin' not in gas:
                self.spectraBottleD[gas], self.spectraBottleD['window_'+gas], self.delayD = self.inter_window(self.spectraBottleD[gas], self.delayD, windowSize=windowSize, direction=direction)
            else:
                self.spectraBottleD[gas], self.spectraBottleD['window_rebin_'+gas], self.rebin_delay = self.inter_window(self.spectraBottleD[gas], self.rebin_delay, windowSize=windowSize, direction=direction)
            
    def inter_window(self, data, delay, windowSize=0, direction='left'):
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
            #self.paddedDelay = np.arange(self.delayB[0])
        # if (stepNum_left%2) == 0:
        #    paddingWindow = np.concatenate(( np.zeros(int(stepNum_left/2))+1,
        #    np.zeros(windowSize*2),np.zeros(int(stepNum_left/2))+1,np.zeros(paddingSize)))
        #    #self.paddedDelay = np.arange(self.delayB[0])
        # else:
        #    paddingWindow = np.concatenate((np.zeros(paddingSize), np.zeros(int((stepNum_left-1)/2))+1,
        #    np.zeros(windowSize*2),np.zeros(int((stepNum_left+1)/2))+1, np.zeros(paddingSize)))
            # plt.plot(dcShift)
        #window=polynomial(window, order=3, plot=False)
        hwindow = np.hanning(len(window))#sps.windows.flattop(len(window),sym=FALSE)
        window2=window#*hwindow #https://stackoverflow.com/questions/55654699/how-to-get-correct-phase-values-using-numpy-fft
        
        return window, window2, delay

    def padding(self, paddingSize=0):
        for gas in list(self.spectraBottleB):
            if "rebin" not in gas:
                self.spectraBottleB[gas], self.delayB = self.inter_padding(self.spectraBottleB[gas], self.delayB, paddingSize = paddingSize)
            else:
                self.spectraBottleB[gas], self.rebin_delay = self.inter_padding(self.spectraBottleB[gas], self.rebin_delay, paddingSize = paddingSize)
        for gas in list(self.spectraBottleD):
            if "rebin" not in gas:
                self.spectraBottleD[gas], self.delayD = self.inter_padding(self.spectraBottleD[gas], self.delayD, paddingSize = paddingSize)
            else:
                self.spectraBottleD[gas], self.rebin_delay = self.inter_padding(self.spectraBottleB[gas], self.rebin_delay, paddingSize = paddingSize)

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
        for gas in self.spectraBottleB.keys():
            iS = np.interp(interpDelay, self.delayB, self.spectraBottleB[gas])
            self.spectraBottleB[gas] = iS
            # print(np.shape(iS))
        self.delayB = interpDelay

    def FFTS(self):
        for gas in self.spectraBottleB.keys():
            if "rebin" not in gas:
                delay = self.delayB
                self.fftSB[gas] = self.FFT(
                    delay, self.spectraBottleB[gas])
            else:
                delay = self.rebin_delay
                self.fftSB[gas] = self.FFT(
                    delay, self.spectraBottleB[gas])
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
        for gas in self.spectraBottleB.keys():
            y = self.spectraBottleB[gas]
            a = np.array([[self.delayB[0], 1], [self.delayB[-1], 1]])
            b = np.array([y[0],y[-1]])
            [k, b] = np.linalg.solve(a, b)
            self.spectraBottleB[gas] = self.spectraBottleB[gas]-(k*self.delayB+b)
            #self.spectraBottleB[gas]=polynomial(self.spectraBottleB[gas], order=3, plot=False)
            hwindow = np.hanning(len(y))#sps.windows.flattop(len(window),sym=FALSE)
            self.spectraBottleB[gas]=self.spectraBottleB[gas]*hwindow

    def runingAverage(self, n=5):
        def runAve(x, N): return np.convolve(x, np.ones(N)/N, mode='valid')
        for gas in self.spectraBottleB.keys():
            self.spectraBottleB[gas] = runAve(self.spectraBottleB[gas], n)
            new_size = len(self.spectraBottleB[gas])
        self.delayB = self.delayB[:new_size]

    def smooth(self, windowSize=100, order=3):
        for gas in self.spectraBottleB.keys():    
            self.spectraBottleB[gas]=sps.savgol_filter(self.spectraBottleB[gas], windowSize, order) # window size 51, polynomial order 3
    
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
        for gas in list(self.spectraBottleB.keys()):
            self.spectraBottleB['rebin_'+gas] = self.rebin_factor(self.spectraBottleB[gas], factor)
        self.rebin_delay = self.rebin_factor(self.delayB, factor)

    def dataProcess(self):
        #self.mass_spectra()
        self.fftSB()

    def useFilter(self, lowcut, highcut):
        for gas in list(self.spectraBottleB.keys()):
            # if 'rebin' not in gas:
            fs = 1/(self.delayB[1]-self.delayB[0])
            # else:
            #     fs = 1/(self.rebin_delay[90]-self.rebin_delay[89])
            self.spectraBottleB['filter_'+gas] = self.butter_bandpass_filter(self.spectraBottleB[gas], lowcut, highcut, fs)

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

    def calDrift(self, _cRange, gas='Ch8'):
        _iS = np.array(np.zeros(self.interSpectraBottleB[gas].shape))
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _iS[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], (3645.84333-1)/33.35641*1e12, (3645.84333+1)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
        _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=2000)
        indexMax = []
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange)
            plt.plot(_delayRange,_inter)
            indexMax = indexMax + [sps.argrelextrema(_inter, np.greater)[0][-1]]
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        #plt.show()

        _ref = sum(indexMax[int(self.interSpectraBottleB[gas].shape[0]/2)-5:int(self.interSpectraBottleB[gas].shape[0]/2)+5])/10
        _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/2000
        #for i in range(self.interSpectraBottleB[gas].shape[0]):
        #    _inter=iS.ev(i,_delayRange+_shift[i])
        #    plt.plot(_delayRange,_inter)
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        #plt.show()
        return _shift

    def delayCorrection(self, _cRange = [400,400+4]):
        for k in range(self.num):
            for gas in self.interSpectraBottleB.keys():
                self.interSpectraBottleB[gas] = self.interinterSpectraBottleB[gas][str(k)]
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

            _delayRange = np.linspace(start=0,stop=1300, num=13000)

            for gas in self.interSpectraBottleB.keys():
                iS = sp.interpolate.RectBivariateSpline(range(self.interSpectraBottleB[gas].shape[0]), self.delayB*1e15, self.interSpectraBottleB[gas])
                inter = 0
                for i in range(self.interSpectraBottleB[gas].shape[0]):
                    inter = inter + iS.ev(i, _delayRange+_shift[i])
                self.spectraBottleB[gas][k] = inter
        self.delayB = _delayRange * 1e-15
        self.interSpectraBottleB = self.spectraBottleB.copy()

        #if self.interSpectraBottleB['Ch0'].shape[0]>1:  
        #    while True:
        #        try:
        #            _shift = self.calDrift(_cRange = _cRange)
        #            break
        #        except IndexError:
        #            print('Wrong Region! Try again!')
        #            _cRange = np.array(_cRange)+0.5
        #    for gas in self.interSpectraBottleB.keys():
        #        iS = sp.interpolate.RectBivariateSpline(range(self.interSpectraBottleB[gas].shape[0]), self.delayB*1e15, self.interSpectraBottleB[gas])
        #        inter = 0
        #        for i in range(self.interSpectraBottleB[gas].shape[0]):
        #            inter = inter + iS.ev(i, _delayRange+_shift[i])
        #        self.interSpectraBottleB[gas]=inter
        #else:
        #    for gas in self.interSpectraBottleB.keys():
        #        self.interSpectraBottleB[gas]=self.interSpectraBottleB[gas].flatten()

    def findZeroDelay(self):
        zeroIndex = {}
        zeroIndex2=[]
        for gas in list(self.interSpectraBottleB.keys()):
            if gas not in ['Ch0', 'Ch2', 'Ch4', 'Ch6', 'Ch8', 'Ch10']:
                continue
            _Spec = self.interSpectraBottleB[gas]
            _Spec=self.butter_bandpass_filter(_Spec, 1/33.35641*1e12, 150/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            plt.plot(_Spec)
            plt.show()
            if gas in ['Ch0', 'Ch2', 'Ch4', 'Ch6']:
                zeroIndex[gas] =(-1*_Spec[:500]).argmax()
                zeroIndex2 = zeroIndex2 + [(-1*_Spec[:500]).argmax()]
            else:
                zeroIndex[gas] =_Spec[:500].argmax()
                zeroIndex2 = zeroIndex2 + [_Spec[:500].argmax()]
        return zeroIndex, zeroIndex2

    def show_FFT(self):
        self.dcRange=0
        gs = gridspec.GridSpec(2, 3)
        fig = plt.figure(figsize=(15,5))
        lab=['Mass 1','Mass 2','Mass 16','Mass 17','Mass 18','ref2']
        #threshold=[0,0,0,0.1,0]#one hour measurement
        #height=[1,0.5,0.5,3,30]
        #threshold=[0.1,0,0,0,0]
        #height=[10,10,10,10,10]
        #distance=[5,5,5,5,30]
        i=0
        if ifsave:
            self.wks = op.new_sheet('w',lname=str('FFT')+str('_')+self.molecule+str('_')+str('%.1E' % Decimal(self.intensity))+str('_')+self.saveRef)
        for gas in self.fftSB.keys():
            if 'filter' in gas  or 'window' in gas or 'rebin' in gas or 'com' in gas:
                continue
            if self.ChN == 11:
                if gas not in ['Ch0','Ch2','Ch4','Ch6','Ch8','Ch10']:
                    continue
                #if gas not in ['Ch1','Ch3','Ch5','Ch7','Ch9']:
                #continue

            else:
                pass

            ax = fig.add_subplot(gs[i])
            axP = ax.twinx()
            label=lab[i]
            #[f, Y, _] = self.fftSB[gas]
            #[f_filter, Y_filter, P_filter] = self.fftSB['filter_'+gas]
            [f_window, Y_window, P_window] = self.fftSB['window_'+gas]
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
            P_window = P_window[bb-aa:bb]

            P_window = np.sum([np.where(np.logical_and(f_window>f-20,f_window<f+20),P_window,0)
             for f in [525.62772,805.21693,883.50191,1338.30036,1591.79458,1696.17455,2147.24515,3205.9563,3645.84333,2669.14501]],0)
            #P_window = np.where(Y_window>0.2,P_window,0)
            #Y_filter=np.where(np.logical_and(f>=10, f<=4500),Y_filter,0)/plotRatio
            #ax.plot(f,np.abs(Y), label=label)#gas)#, label=self.trueScanTime)
            #ax.plot(f,np.abs(Y_filter)/20, label='filter'+label)#gas)#, label=self.trueScanTime)
            ax.plot(f_window, np.abs(Y_window), label=label)#gas)#, label=self.trueScanTime)
            #ax.plot(f_rebin,np.abs(Y_rebin), label='rebin_window_'+label)#gas)#, label=self.trueScanTime)
            #ax.plot(f,np.real(Y), label="Re")#gas)#, label=self.trueScanTime)
            #ax.plot(f,np.imag(Y), label="Im")#gas)#, label=self.trueScanTime)
            axP.plot(f_window,P_window/np.pi,'r')

            #ax.plot(f[peaksIndex[i]],np.abs(Y_window)[peaksIndex[i]], "rx")
            #if gas == 'Ch2+4':
            #    [f_window, Y_ref, P_window] = self.fftSB['window_'+'Ch2-4']
            #    ax.plot(f, np.abs(Y_ref), label='Ch2-4')
            ax.set_xlabel('Frequency/cm-1')
            ax.set_ylabel('Amplitude')
            #axP.set_ylabel('Phase/pi')
            ax.legend(loc=2)
            #ax.set_ylim([0,15])

            ax.set_xlim([450,4500])
            #ax.set_xlim([58870,59920])
            #ax.set_title(gas)
            i=i+1
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('FFT')+label+self.molecule+str('%.1E' % Decimal(self.intensity))+str('.dat')), np.abs(Y_window))
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('FFT')+label+self.molecule+str('%.1E' % Decimal(self.intensity))+str('frequency')+str('.dat')), f)
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
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('FFT')+str(gas)+self.molecule+str('%.1E' % Decimal(self.intensity))+str('.dat')), np.abs(Y))
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('FFT')+str(gas)+self.molecule+str('%.1E' % Decimal(self.intensity))+str('frequency')+str('.dat')), f*33.35641)
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
        lab=['Mass1','Mass2','Mass16','Mass17','Mass18','Ch2+Ch4','Mass1+Mass16', 'Mass1+Mass17', 'Mass16+Mass17','Mass17+Mass18']
        i=0
        if ifsaveT:
            self.wks = op.new_sheet('w',lname=str('Time')+str('_')+self.molecule+str('_')+str('%.1E' % Decimal(self.intensity))+str('_')+self.saveRef)
        for gas in self.spectraBottleB.keys():
            print(gas)
            #if 'rebin' not in gas or 'filter' not in gas:
            if 'filter' in gas or 'window' in gas or 'rebin' in gas or 'com' in gas:
                continue
            #elif 'rebin' not in gas:
            
            if self.ChN == 11:
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
                     self.spectraBottleB["window_"+gas], label=gas)#/np.amax(self.spectraBottleB[gas])
            #ax.plot(delay*10**15+shift,
            #         self.spectraBottleB['window_filter_'+gas], label=label)
            #ax.plot(delay*10**15+shift,
            #         self.spectraBottleB['window_'+gas], label='window_'+gas)
            #ax.plot(rebindelay*10**15+shift,
            #         self.spectraBottleB['rebin_window_'+gas], label='rebin_window_'+gas)
            ax.set_xlabel("Delay/fs")
            ax.set_ylabel('a.u.')
            #plt.xlim([300,400])
            ax.legend()
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('Time')+label+self.molecule+self.stage+str('%.1E' % Decimal(self.intensity))+str('.dat')), x.spectraBottleB[gas])
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('delay')+str('%.1E' % Decimal(self.intensity))+str('.dat')), delay*10**15)
            if ifsaveT:
                self.wks.from_list(0, delay, 'X')
                self.wks.from_list(i, self.spectraBottleB[gas], lname=gas, axis='Y')
        #plt.legend()
        fig.tight_layout()
        #plt.savefig("spectra_31"+d,dpi=720, bbox_inches='tight', pad_inches=0.2)
        plt.show()


    def phaseRetrive2(self,omega=[3645.84333]):
        self.filterSpecMass1 = {}
        self.filterSpecMass2 = {}
        self.filterSpecMass16 = {}
        self.filterSpecMass17 = {}
        self.filterSpecMass18 = {}
        self.filterPhaseMass1 = {}
        self.filterPhaseMass2 = {}
        self.filterPhaseMass16 = {}
        self.filterPhaseMass17 = {}
        self.filterPhaseMass18 = {}
        for f in omega:
            f0=(f-0.1)/33.35641*1e12
            f1=(f+0.1)/33.35641*1e12
            self.filterSpecMass1[str(f)] = self.useFilter2(self.spectraBottleB_noPad["Ch0"],f0,f1)
            self.filterSpecMass2[str(f)] = self.useFilter2(self.spectraBottleB_noPad["Ch2"],f0,f1)
            self.filterSpecMass16[str(f)] = self.useFilter2(self.spectraBottleB_noPad["Ch4"],f0,f1)
            self.filterSpecMass17[str(f)] = self.useFilter2(self.spectraBottleB_noPad["Ch6"],f0,f1)
            self.filterSpecMass18[str(f)] = self.useFilter2(self.spectraBottleB_noPad["Ch8"],f0,f1)
            self.filterPhaseMass1[str(f)] = np.angle(sps.hilbert(self.filterSpecMass1[str(f)]))
            self.filterPhaseMass2[str(f)] = np.angle(sps.hilbert(self.filterSpecMass2[str(f)]))
            self.filterPhaseMass16[str(f)] = np.angle(sps.hilbert(self.filterSpecMass16[str(f)]))
            self.filterPhaseMass17[str(f)] = np.angle(sps.hilbert(self.filterSpecMass17[str(f)]))
            self.filterPhaseMass18[str(f)] = np.angle(sps.hilbert(self.filterSpecMass18[str(f)]))
            sampleIndex = np.zeros(self.delayB_noPad.size)+1
            #print(str(f)+": "+str(self.filterPhaseMass1[str(f)][1000]/np.pi)+str(", ")+str(self.filterPhaseMass2[str(f)][1000]/np.pi)+str(", ")+str(self.filterPhaseMass16[str(f)][1000]/np.pi)+str(", ")+str(self.filterPhaseMass17[str(f)][1000]/np.pi)+str(", ")+str(self.filterPhaseMass18[str(f)][1000]/np.pi))
            sampleIndex = np.where(np.abs(self.filterPhaseMass18[str(f)]-np.pi)<2*np.pi*2/9.15*(self.delayB[1]-self.delayB[0])*1e15,sampleIndex,0)
            sampleIndex = np.where(self.delayB_noPad*1e15<800, sampleIndex, 0)
            sampleIndex = np.where(self.delayB_noPad*1e15>200, sampleIndex, 0)
        
            #plt.plot(sampleIndex)
            #plt.show()
            #plt.plot(self.filterSpecMass1[str(f)])
            #plt.plot(self.filterSpecMass17[str(f)])
            #plt.plot(self.filterSpecMass18[str(f)])
            #plt.plot(self.filterPhaseMass1[str(f)],label='Mass1')
            #plt.plot(self.filterPhaseMass17[str(f)],label='Mass17')
            #plt.plot(self.filterPhaseMass18[str(f)],label='Mass18')
            #plt.legend()
            #plt.xlabel('delay')
            #plt.ylabel('Phase')
            #plt.show()
            self.retrivedPhase[str(f)] = [np.sum(self.filterPhaseMass1[str(f)]*sampleIndex), np.sum(self.filterPhaseMass2[str(f)]*sampleIndex), np.sum(self.filterPhaseMass16[str(f)]*sampleIndex),
             np.sum(self.filterPhaseMass17[str(f)]*sampleIndex), np.sum(self.filterPhaseMass18[str(f)]*sampleIndex)]/np.sum(sampleIndex)
        return 

if __name__ == '__main__':
    
    f0 = [r'scan_tof_2022-05-24-17-38-18']
    f1=[r'scan_tof_2022-05-25-09-47-24']
    f2=[r'scan_tof_2022-05-25-12-31-23']
    f3=[r'scan_tof_2022-05-25-14-08-02']
    f4=[r'scan_tof_2022-05-25-15-50-01']
    f5=[r'scan_tof_2022-05-25-17-27-19']
    f6=[r'scan_tof_2022-05-25-18-56-51']
    f7=[r'scan_tof_2022-05-26-12-44-27']
    f8=[r'scan_tof_2022-05-26-14-13-58']
    f9=[r'scan_tof_2022-05-26-15-46-15']
    f10=[r'scan_tof_2022-05-26-17-00-01',r'scan_tof_2022-05-26-22-52-13']
    f11=[r'scan_tof_2022-05-26-18-50-20']
    f12=[]
    f13=[r'scan_tof_2022-05-27-10-23-47']
    f14=[r'scan_tof_2022-05-27-11-54-12',r'scan_tof_2022-05-27-19-53-05']
    f15=[r'scan_tof_2022-05-27-13-12-13']
    f16=[r'scan_tof_2022-05-27-15-11-46']
    f17=[r'scan_tof_2022-05-27-16-32-58']
    f18=[r'scan_tof_2022-05-27-18-03-09']

    #since 20220530
    
    f19 = [r'scan_tof_2022-05-30-10-43-15']
    f20 = [r'scan_tof_2022-05-30-12-23-32']
    f21 = [r'scan_tof_2022-05-30-13-50-56']
    f22 = [r'scan_tof_2022-05-30-15-27-08']
    f23 = [r'scan_tof_2022-05-30-16-57-23']
    f24 = [r'scan_tof_2022-05-30-18-18-33']
    f25 = [r'scan_tof_2022-05-30-19-44-25']
    f26 = [r'scan_tof_2022-05-31-09-23-21']
    f27 = [r'scan_tof_2022-05-31-11-12-23']
    f28 = [r'scan_tof_2022-05-31-12-26-11']    
    f29 = [r'scan_tof_2022-05-31-14-16-24']
    f30 = [r'scan_tof_2022-05-31-15-42-05']
    f31 = [r'scan_tof_2022-05-31-17-51-43']
    f32 = [r'scan_tof_2022-05-31-19-41-25']
    #1st June, today the fiber output is a little bit less, and the data also shows lower intensity though try to set as same parameter as last round.
    f33 = [r'scan_tof_2022-06-01-10-40-05']
    f34 = [r'scan_tof_2022-06-01-12-10-05']
    f35 = [r'scan_tof_2022-06-01-13-39-07']
    f36 = [r'scan_tof_2022-06-01-15-20-40']
    f37 = [r'scan_tof_2022-06-01-16-54-40']
    #r'scan_tof_2022-06-04-22-37-09',r'scan_tof_2022-06-01-18-10-15',r'scan_tof_2022-06-04-12-44-47',r'scan_tof_2022-06-01-15-20-40',r'scan_tof_2022-06-01-20-53-44',r'scan_tof_2022-06-02-17-31-08',
    f38 = [r'scan_tof_2022-06-04-22-37-09', r'scan_tof_2022-06-04-15-39-25']#r'scan_tof_2022-06-04-22-37-09', r'scan_tof_2022-06-04-15-39-25'
    f38r =[r'scan_tof_2022-06-03-22-44-54',r'scan_tof_2022-06-04-18-15-48',r'scan_tof_2022-06-01-20-53-44',r'scan_tof_2022-06-04-12-44-47']
    f38rr=[r'scan_tof_2022-06-04-15-39-25',r'scan_tof_2022-06-03-22-44-54',r'scan_tof_2022-06-04-18-15-48',r'scan_tof_2022-06-04-12-44-47']
    f38rrr = [r'scan_tof_2022-06-04-22-37-09',r'scan_tof_2022-06-01-18-10-15',r'scan_tof_2022-06-02-18-58-56',r'scan_tof_2022-06-03-09-31-11',r'scan_tof_2022-06-03-12-20-18',r'scan_tof_2022-06-03-18-10-43',r'scan_tof_2022-06-01-15-20-40',r'scan_tof_2022-06-02-17-31-08']
    #r'scan_tof_2022-06-04-22-37-09',r'scan_tof_2022-06-01-18-10-15',r'scan_tof_2022-06-02-18-58-56',r'scan_tof_2022-06-03-09-31-11',,r'scan_tof_2022-06-03-12-20-18',r'scan_tof_2022-06-03-18-10-43',r'scan_tof_2022-06-01-15-20-40',r'scan_tof_2022-06-02-17-31-08'
    f39 = [r'scan_tof_2022-06-02-10-12-43']
    f43 = [r'scan_tof_2022-06-02-17-31-08']
    f44 = [r'scan_tof_2022-06-02-18-58-56']
    f45 = [r'scan_tof_2022-06-03-09-31-11',r'scan_tof_2022-06-03-10-56-37']
    f46 = [r'scan_tof_2022-06-03-16-24-18']

    f47 = [r'scan_tof_2022-06-05-15-07-15',r'scan_tof_2022-06-05-13-10-09']
    f48 = [r'scan_tof_2022-06-06-19-53-17',r'scan_tof_2022-06-06-15-20-36',r'scan_tof_2022-06-05-20-18-21',r'scan_tof_2022-06-06-13-20-33',] #200-200
    #r'scan_tof_2022-06-06-18-27-03',r'scan_tof_2022-06-05-17-03-06', 200-200
    #r'scan_tof_2022-06-07-09-56-51' 200-200 delay drift ...try to relocate the zero delay later
    #scan_tof_2022-06-07-11-29-36 200-200 delay little drift ...try to relocate the zero delay later
    #scan_tof_2022-06-07-13-33-30 250-230
    #f49 = [r'scan_tof_2022-06-07-16-29-08']
    f49 = [r'scan_tof_2022-06-07-17-41-15']
    ftest = [r'scan_tof_2022-06-22-19-24-23']
    #250-230 good scan_tof_2022-06-07-13-33-30 scan_tof_2022-06-07-16-29-08
    ##250-230 bad scan_tof_2022-06-07-17-41-15

    #phaseRetrive
    fp1=[r'scan_tof_2022-06-04-22-37-09']
    fp2=[r'scan_tof_2022-06-04-15-39-25',r'scan_tof_2022-06-03-22-44-54',]
    fp3=[r'scan_tof_2022-06-04-18-15-48']
    fp4 = fp1+fp2+fp3


    ############################20220621
    f50 = [r'scan_tof_2022-06-23-10-09-03',r'scan_tof_2022-06-23-11-48-23']
    f51 = [r'scan_tof_2022-06-23-13-11-05'] #190-200
    f52 = [r'scan_tof_2022-06-23-15-25-48'] #180-200
    f53 = [r'scan_tof_2022-06-23-17-22-11',r'scan_tof_2022-06-24-11-11-20',r'scan_tof_2022-06-23-20-24-23',r'scan_tof_2022-06-24-13-33-25 - Copy',r'scan_tof_2022-06-24-15-11-10',r'scan_tof_2022-06-24-16-40-40',r'scan_tof_2022-06-24-20-31-22'] #180-180
    #f53 = [r'scan_tof_2022-06-23-20-24-23',r'scan_tof_2022-06-24-13-33-25 - Copy'] #180-180
    #f53 = [r'scan_tof_2022-06-24-11-11-20 - Copy',r'scan_tof_2022-06-24-13-33-25 - Copy'] #180-180
    #f53 = [r'scan_tof_2022-06-24-15-11-10'] #180-180
    #f53 = [r'scan_tof_2022-06-24-20-31-22']
    #f53 = [r'scan_tof_2022-06-27-11-08-21 - Copy']#zero delay drift
    #f53 = [r'scan_tof_2022-06-27-13-40-24']#
    #f53 = [r'scan_tof_2022-06-27-17-09-33']
    f54 = [r'scan_tof_2022-06-27-19-55-38']
    f54 = [r'scan_tof_2022-06-28-12-45-25']
    f54 = [r'scan_tof_2022-06-28-15-47-53']
    #f54 = [r'scan_tof_2022-06-28-19-18-27']
    #f54 = [r'scan_tof_2022-06-29-12-09-27']
    f55 = [r'scan_tof_2022-06-29-18-16-59']
    f56 = [r'scan_tof_2022-07-03-12-09-43']
    
    
    
    
    
    
   


    fp1=FFT_ionS(filename=fp1, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    fp2=FFT_ionS(filename=fp2, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    fp3=FFT_ionS(filename=fp3, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    fp4=FFT_ionS(filename=fp4, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)

    fall = f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18+f19+f20+f21+f22+f23+f24+f25+f26+f27+f28+f29+f30+f31+f32+f33+f34+f35+f36+f37+f38+f39+f43+f44+f45+f46+ftest
    
    #probe 230 mW
    f0 = FFT_ionS(filename=f0, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)
    f1 = FFT_ionS(filename=f1, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)
    f2 = FFT_ionS(filename=f2, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(60,8.2,5),ChN=11)
    f3 = FFT_ionS(filename=f3, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11)
    f4 = FFT_ionS(filename=f4, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)
    f5 = FFT_ionS(filename=f5, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)#repeat

    f19 = FFT_ionS(filename=f19, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)
    f20 = FFT_ionS(filename=f20, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.5,5),ChN=11)#maybe repeat again
    f21 = FFT_ionS(filename=f21, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.5,5),ChN=11)
    f22 = FFT_ionS(filename=f22, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.5,5),ChN=11)

    #probe 100 mW
    f6 = FFT_ionS(filename=f6, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)
    f16 = FFT_ionS(filename=f16, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)
    f17 = FFT_ionS(filename=f17, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)#repeat
    f18 = FFT_ionS(filename=f18, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11)


    #probe 150 mW
    f7 = FFT_ionS(filename=f7, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)
    f13 = FFT_ionS(filename=f13, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)
    f14 = FFT_ionS(filename=f14, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.5,5),ChN=11)#repeat
    f15 = FFT_ionS(filename=f15, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11)
    
    f34 = FFT_ionS(filename=f34, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)
    f35 = FFT_ionS(filename=f35, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9.8,5),ChN=11)#repeat 20220601 actual intensity seems lower due to longer pulse
    f36 = FFT_ionS(filename=f36, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)#good repeat 20220601
    f37 = FFT_ionS(filename=f37, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,11.6,5),ChN=11)
    #f38 = FFT_ionS(filename=f38, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    f45 = FFT_ionS(filename=f45, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11) #actua power 155-155
    
    #probe 180mW
    f8 = FFT_ionS(filename=f8, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)#repeat
    f9 = FFT_ionS(filename=f9, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)#interesting, multi fre in H2+ 
    f10 = FFT_ionS(filename=f10, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)#Mass2 and Mass 16 exactly anti phase
    f11 = FFT_ionS(filename=f11, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11)
    #f12 = FFT_ionS(filename=f12, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)
    
    f29 = FFT_ionS(filename=f29, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11) #MCP 3600
    f30 = FFT_ionS(filename=f30, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11) #MCP 3800
    f31 = FFT_ionS(filename=f31, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)
    f32 = FFT_ionS(filename=f32, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11) #3-hour measurement
    f33 = FFT_ionS(filename=f33, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11) #repeat 20220601

    #probe 200mW
    f23 = FFT_ionS(filename=f23, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)
    f24 = FFT_ionS(filename=f24, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.5,5),ChN=11)
    f25 = FFT_ionS(filename=f25, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.5,5),ChN=11) #MCP 3200V
    f26 = FFT_ionS(filename=f26, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.5,5),ChN=11) #MCP 3500V 1.3e-7
    f27 = FFT_ionS(filename=f27, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.5,5),ChN=11)
    f28 = FFT_ionS(filename=f28, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.5,5),ChN=11)
    fall = FFT_ionS(filename=fall, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)

    #same pump probe power
    f38 = FFT_ionS(filename=f38, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    f38r = FFT_ionS(filename=f38r, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    f38rr = FFT_ionS(filename=f38rrr, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    f39 = FFT_ionS(filename=f39, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(170,10.7,5),ChN=11)#laser drift
    ftest = FFT_ionS(filename=ftest, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)

    #try
    f43 = FFT_ionS(filename=f43, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(170,10.7,5),ChN=11)#150-140
    f44 = FFT_ionS(filename=f44, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(170,10.7,5),ChN=11)#150-160
    f46 = FFT_ionS(filename=f46, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(170,10.7,5),ChN=11)#150-160
    f47 = FFT_ionS(filename=f47, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(180,10.7,5),ChN=11)#180-180
    f48 = FFT_ionS(filename=f48, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9.5,5),ChN=11)#200-200
    f49 = FFT_ionS(filename=f49, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)#200-200

    f50 = FFT_ionS(filename=f50, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9,5),ChN=11)#200-200
    f51 = FFT_ionS(filename=f51, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(190,9,5),ChN=11)#190-200
    f52 = FFT_ionS(filename=f52, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(180,9,5),ChN=11)#180-200
    f53 = FFT_ionS(filename=f53, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9,5),ChN=11)#180-180
    f54 = FFT_ionS(filename=f54, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9,5),ChN=11)#200-200
    f55 = FFT_ionS(filename=f55, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9,5),ChN=11)#150-150
    f56 = FFT_ionS(filename=f56, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9,5),ChN=11)#200-150


    #ff=[fp1,fp2,fp3,fp4]
    ff=[f48]
    for x in ff:
        x.pathFinder()
        print(x.filepath)
        
        x.read_splitB(isCorrection=True,sumN=1)
        #x.findZeroDelay()
        #x.window(windowSize=150, direction='left')
        #x.rmvExp()
        #x.smooth(windowSize=9)
        #x.show_Spectra()
        #plt.show()
        #x.padding(paddingSize=100000)
        #x.FFTS()
        #x.show_FFT()
        #x.phaseRetrive2()
        
    