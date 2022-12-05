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
    def __init__(self, filename, scanTime, sampRate, massRange=[5, 50], molecule='H2', intensity = 0, dcRange=5, ChN=10, cal_mass=[2,16], cal_pixel=[3082,10216]):
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
        self.fftSB = {}
        self.stftSB = {}
        self.waveletsCoe = {}
        self.zeroIndex = {}
        self.zeroIndexMin = 10000000000000000000000000
        self.zeroIndexMax = 0
        self.spec = None
        self.scanNo = None
        self.ChN = ChN
        self.molecule = molecule
        self.intensity = intensity
        self.stage = 'piezo'

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
            self.exportPath = pl.PureWindowsPath(self.rootPath, r'plot20220602')
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

    def read_splitB(self, isCorrection = False, sumN=10):
        ChN=self.ChN #number of channels measured
        for fp in self.filepath:
            if isCorrection:
                self.delayCorrection(n=sumN, fp=fp)#n means correct every n scan
            else:

                for i in [0,2,4,6,8,10]:
                    self.spectraBottleB['Ch'+str(i)] = 0
            
                with h5py.File(fp, 'r+') as f:
                    print(f.keys())
                    data = np.array(f['dataB'], dtype=float)
                    print(data.shape)
                    for i in [0,2,4,6,8,10]:
                            self.interSpectraBottleB['Ch'+str(i)] = np.sum(data[i::ChN][:350],0)
                            #iii=data[i::ChN]
                            #for i in range(1,5):
                            #    plt.plot(iii[i*3])
                            #plt.show()
                            #self.spectraBottleB['Ch'+str(i)] = self.interSpectraBottleB['Ch'+str(i)]+self.spectraBottleB['Ch'+str(i)]

            self.zeroIndex[str(fp)], __inter = self.findZeroDelay()
            self.zeroIndexMin = np.amin(__inter+[self.zeroIndexMin])
            self.zeroIndexMax = np.amax(__inter+[self.zeroIndexMax])

        for fp in self.filepath:
            if isCorrection:
                self.delayCorrection(n=sumN, fp=fp)#n means correct every n scan
                for i in [0,2,4,6,8,10]:
                    self.spectraBottleB['Ch'+str(i)] = self.spectraBottleB['Ch'+str(i)] + self.interSpectraBottleB['Ch'+str(i)][(self.zeroIndex[str(fp)]['Ch'+str(i)]-self.zeroIndexMin):(self.scanLengthB-(self.zeroIndexMax-self.zeroIndex[str(fp)]['Ch'+str(i)]))]
                    #plt.plot( self.spectraBottleB['Ch'+str(i)] )
                    #plt.show()
            else:
                with h5py.File(fp, 'r+') as f:
                    print(f.keys())
                    data = np.array(f['dataB'], dtype=float)
                    print(data.shape)
                    #for i in range(ChN):
                    for i in [0,2,4,6,8,10]:
                            self.interSpectraBottleB['Ch'+str(i)] = np.sum(data[i::ChN],0)
                            #self.interSpectraBottleB['Ch'+str(i)] = self.interSpectraBottleB['Ch'+str(i)]#/np.mean(self.interSpectraBottleB['Ch'+str(i)][1000:3000])
                            self.spectraBottleB['Ch'+str(i)] = self.spectraBottleB['Ch'+str(i)] + self.interSpectraBottleB['Ch'+str(i)][(self.zeroIndex[str(fp)]['Ch'+str(i)]-self.zeroIndexMin):(self.scanLengthB-(self.zeroIndexMax-self.zeroIndex[str(fp)]['Ch'+str(i)]))]
                            #plt.plot( self.spectraBottleB['Ch'+str(i)] )
                            #plt.show()
                            #fix the length of the trace
                            #if fn == 0:
                            #    self.spectraBottleB['Ch'+str(i)] = 1*self.spectraBottleB['Ch'+str(i)]#/np.mean(self.spectraBottleB['Ch'+str(i)][1000:3000])
                    #for i in range(ChN):
                    #    if np.sum(self.spectraBottleB['Ch'+str(i)])==0 or self.interSpectraBottleB['Ch'+str(i)].size==self.spectraBottleB['Ch'+str(i)].size:
                    #       self.spectraBottleB['Ch'+str(i)] = self.spectraBottleB['Ch'+str(i)] + self.interSpectraBottleB['Ch'+str(i)]
                    #    elif self.interSpectraBottleB['Ch'+str(i)].size > self.spectraBottleB['Ch'+str(i)].size:
                    #        interL=self.spectraBottleB['Ch'+str(i)].size
                    #        self.spectraBottleB['Ch'+str(i)] = self.spectraBottleB['Ch'+str(i)] + self.interSpectraBottleB['Ch'+str(i)][:interL]
                    #    else:
                    #        interL=self.interSpectraBottleB['Ch'+str(i)].size
                    #        self.spectraBottleB['Ch'+str(i)] = self.spectraBottleB['Ch'+str(i)][:interL] + self.interSpectraBottleB['Ch'+str(i)]

        for i in [0,2,4,6,8,10]:
            self.spectraBottleB['Ch'+str(i)] = self.spectraBottleB['Ch'+str(i)][:3000]
        delayStep = self.delayB[1]-self.delayB[0]
        
        self.delayB = (np.arange(self.spectraBottleB['Ch'+str(0)].size)-self.zeroIndexMin)*delayStep
        self.spectraBottleB['Ch'+str(2)] = self.spectraBottleB['Ch'+str(2)]-self.spectraBottleB['Ch'+str(10)]
        self.spectraBottleB['Ch'+str(8)] = self.spectraBottleB['Ch'+str(8)]-self.spectraBottleB['Ch'+str(4)]
        
        add = lambda i,j, m: self.spectraBottleB['Ch'+str(i)]+self.spectraBottleB['Ch'+str(j)]*m


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
        self.spectraBottleB_noPad=self.spectraBottleB.copy()
        self.delayB_noPad = self.delayB.copy()
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

    def interpS(self, interNum=1):
        gas='Ch8'
        interpDelay = np.arange(
            self.delayB[0], self.delayB[-1], (self.delayB[1]-self.delayB[0])/(interNum+1))
        iS = sp.interpolate.RectBivariateSpline(range(self.interSpectraBottleB[gas].shape[0]), self.delayB, self.interSpectraBottleB[gas])
        # print(np.shape(iS))

    def findZeroDelay(self):
        zeroIndex = {}
        zeroIndex2=[]
        for gas in list(self.interSpectraBottleB.keys()):
            if gas not in ['Ch0', 'Ch2', 'Ch4', 'Ch6', 'Ch8', 'Ch10']:
                continue
            _Spec = self.interSpectraBottleB[gas]
            _Spec=self.butter_bandpass_filter(_Spec, 1/33.35641*1e12, 150/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            #plt.plot(_Spec)
            #plt.show()
            if gas in ['Ch0', 'Ch2', 'Ch4', 'Ch6']:
                zeroIndex[gas] =(-1*_Spec[:500]).argmax()
                zeroIndex2 = zeroIndex2 + [(-1*_Spec[:500]).argmax()]
            else:
                zeroIndex[gas] =_Spec[:500].argmax()
                zeroIndex2 = zeroIndex2 + [_Spec[:500].argmax()]
        return zeroIndex, zeroIndex2
    
    def findZeroDelayC(self, m):
        zeroIndex = {}
        zeroIndex2=[]
        _range = [100,250]
        
        for gas in list(self.interSpectraBottleB.keys()):
            if gas not in ['Ch0', 'Ch2', 'Ch4', 'Ch6', 'Ch8', 'Ch10']:
                continue
            _Spec = self.interSpectraBottleB[gas][m]
            _Spec=self.butter_bandpass_filter(_Spec, (3645.84333-1)/33.35641*1e12, (3645.84333+1)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            if gas == 'Ch6':    
                plt.plot(_Spec[1000:1018])
                
            if gas in ['Ch0', 'Ch2', 'Ch4']:
                zeroIndex[gas] =_Spec[_range[0]:_range[1]].argmin()
                zeroIndex2 = zeroIndex2 + [_Spec[_range[0]:_range[1]].argmin()]
            else:
                zeroIndex[gas] =_Spec[_range[0]:_range[1]].argmax()
                zeroIndex2 = zeroIndex2 + [_Spec[_range[0]:_range[1]].argmax()]
        
        return zeroIndex, zeroIndex2

    def findZeroDelayD(self):
        zeroIndex = []
        for fp in self.filepath:
            with h5py.File(fp, 'r+') as f:
                indexRange = [int(self.calculator.cal_pixel(0.5)),int(self.calculator.cal_pixel(16.5))]
                interSpec = np.sum(np.array(f['useful'])[indexRange[0]:indexRange[1]],0)
                interSpec=self.butter_bandpass_filter(interSpec, 100/33.35641*1e12, 1000/33.35641*1e12, 1/(self.delayD[1]-self.delayD[0]))

            #plt.plot(interSpec)
            #plt.show()
            zeroIndex = zeroIndex+[np.abs(interSpec).argmax()]
        print(zeroIndex)
        return zeroIndex

    def delayCorrection(self, n, fp):
        ChN=self.ChN #number of channels measured
        for i in [0,2,4,6,8,10]:
            self.spectraBottleB['Ch'+str(i)] = 0

        with h5py.File(fp, 'r+') as f:
            print(f.keys())
            data = np.array(f['dataB'], dtype=float)
            print(data.shape)
            m = int(data.shape[0]/self.ChN/n)
            for i in [0,2,4,6,8,10]:
                    self.interSpectraBottleB['Ch'+str(i)] = data[i::ChN][:int(m*n)].reshape(m, n, data.shape[1]).sum(axis=1)
                    print(self.interSpectraBottleB['Ch'+str(i)].shape)
        #self.interpS(79)
        for j in range(m):
            self.zeroIndex[str(j)], __inter = self.findZeroDelayC(m=j)
            self.zeroIndexMin = np.amin(__inter+[self.zeroIndexMin])
            self.zeroIndexMax = np.amax(__inter+[self.zeroIndexMax])
        print((self.delayB[1]-self.delayB[0])*1e15)
        plt.show()
        for j in range(m):
            for i in [0,2,4,6,8,10]:
                self.spectraBottleB['Ch'+str(i)] = self.spectraBottleB['Ch'+str(i)] + self.interSpectraBottleB['Ch'+str(i)][j][(self.zeroIndex[str(j)]['Ch'+str(i)]-self.zeroIndexMin):(self.scanLengthB-(self.zeroIndexMax-self.zeroIndex[str(j)]['Ch'+str(i)]))]
        for i in [0,2,4,6,8,10]:
            self.interSpectraBottleB['Ch'+str(i)] = self.spectraBottleB['Ch'+str(i)]
            self.spectraBottleB['Ch'+str(i)] = 0

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
        d='_80pp_240pr'
        #plt.savefig('20220214_'+"fft_"+d,dpi=720, bbox_inches='tight', pad_inches=0.05)
        plt.show()

    def show_FFT2(self):
        self.dcRange=0
        gs = gridspec.GridSpec(2, 3)
        fig = plt.figure(figsize=(15,5))
        lab=['Mass 1+17','Mass 1-17','Mass 17+18','Mass 17-18']
        #threshold=[0,0,0,0.1,0]#one hour measurement
        #height=[1,0.5,0.5,3,30]
        #threshold=[0.1,0,0,0,0]
        #height=[10,10,10,10,10]
        #distance=[5,5,5,5,30]
        i=11
        for gas in self.fftSB.keys():
            if 'filter' in gas  or 'window' in gas or 'rebin' in gas or 'com' in gas:
                continue
            if self.ChN == 11:
                if gas not in ['Ch0+6','Ch0-6','Ch6+8','Ch6-8']:
                    continue
                #if gas not in ['Ch1','Ch3','Ch5','Ch7','Ch9']:
                #continue

            else:
                pass

            ax = fig.add_subplot(gs[i])
            #axP = ax.twinx()
            label=lab[i]
            [f, Y, _] = self.fftSB[gas]
            #[f_filter, Y_filter, P_filter] = self.fftSB['filter_'+gas]
            [f_window, Y_window, P_window] = self.fftSB['window_'+gas]
            #[f_rebin, Y_rebin, P_rebin] = self.fftSB['rebin_'+gas]
            plotRatio=1
            f=f.astype(float)*33.35641
            #f_rebin=f_rebin*33.35641
            Y_window=np.where(f>=100,Y_window,0)/np.amax(Y_window)
            #Y_filter=np.where(np.logical_and(f>=10, f<=4500),Y_filter,0)/plotRatio
            #ax.plot(f,np.abs(Y), label=label)#gas)#, label=self.trueScanTime)
            #ax.plot(f,np.abs(Y_filter)/20, label='filter'+label)#gas)#, label=self.trueScanTime)
            ax.plot(f, np.abs(Y_window), label=label)#gas)#, label=self.trueScanTime)
            #ax.plot(f_rebin,np.abs(Y_rebin), label='rebin_window_'+label)#gas)#, label=self.trueScanTime)
            #ax.plot(f,np.real(Y), label="Re")#gas)#, label=self.trueScanTime)
            #ax.plot(f,np.imag(Y), label="Im")#gas)#, label=self.trueScanTime)
            #axP.plot(f,P_filter/np.pi+1,'r')

            #ax.plot(f[peaksIndex[i]],np.abs(Y_window)[peaksIndex[i]], "rx")
            #if gas == 'Ch2+4':
            #    [f_window, Y_ref, P_window] = self.fftSB['window_'+'Ch2-4']
            #    ax.plot(f, np.abs(Y_ref), label='Ch2-4')
            ax.set_xlabel('Frequency/cm-1')
            ax.set_ylabel('Amplitude')
            #axP.set_ylabel('Phase/pi')
            ax.legend(loc=2)
            #ax.set_ylim([0,15])

            ax.set_xlim([0,4500])
            #ax.set_xlim([58870,59920])
            #ax.set_title(gas)
            i=i+1
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('FFT')+label+self.molecule+str('%.1E' % Decimal(self.intensity))+str('.dat')), np.abs(Y_window))
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('FFT')+label+self.molecule+str('%.1E' % Decimal(self.intensity))+str('frequency')+str('.dat')), f)
            if ifsave:    
                self.wks.from_list(0, f, 'X')
                self.wks.from_list(i, np.abs(Y_window), lname=label, axis='Y')

        plt.legend()
        #fig.tight_layout()
        d='_80pp_240pr'
        #plt.savefig('20220214_'+"fft_"+d,dpi=720, bbox_inches='tight', pad_inches=0.05)
        plt.show()

    def show_FFTD(self):
        #plt.figure()
        #unused_f, ref = self.fftS['CH3OH+']
        i=12
        for gas in self.fftSD.keys():
            if 'Ch2' not in gas:
                continue
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

    def showRelativePhase(self, shift=0):
        fig= plt.figure(figsize=(8,3))
        bax = brokenaxes(xlims=((0, 1200),),)#ylims=((-10,5),)
        #lab={'Ch0':'Mass1','Ch1':'Mass2','Ch2':'Mass16 x 10','Ch3':'Mass17 x 0.1','Ch4':'Mass18 x 0.1'}
        lab={'Ch0':'Mass1','Ch1':'Mass2x 4','Ch2':'Mass16 x 8','Ch3':'Mass17 x 8','Ch4':'Mass18 x 15', 'Ch0+1':'Mass1+2'}
        i=0
        for gas in self.spectraBottleB.keys():
            #if 'rebin' not in gas or 'filter' not in gas:
            if 'fil' not in gas or 'window' not in gas or 'rebin' in gas or 'com' in gas:
                continue
            #elif 'Ch1' not in gas:
            #    continue
            label=lab[gas[14:]]
            if 'com' in gas:
                label = 'filter1543-1623cm_'+label
            else:
                label = 'filter3500-4250cm_'+label

            i=i+1
            if 'com' in gas:
                delay = self.delayB*1816/3647
            else:
                delay = self.delayB
            rebindelay = self.rebin_delay
            #else:
            #    delay = self.rebin_delay
            if 'Ch2' in gas or 'Ch3' in gas:
                bax.plot(delay*10**15+shift,
                        self.spectraBottleB[gas]*8,label=label)
            elif 'Ch1' in gas:
                bax.plot(delay*10**15+shift,
                        self.spectraBottleB[gas]*4,label=label)
            elif 'Ch4' in gas:
                bax.plot(delay*10**15+shift,
                        self.spectraBottleB[gas]*15,label=label)
            else:
                bax.plot(delay*10**15+shift,
                        self.spectraBottleB[gas],label=label)
            bax.set_xlabel("Delay/fs")
            bax.set_ylabel('a.u.')
        bax.legend(ncol=2, loc=3)
        #plt.ylim([-1.5,1.5])
        #fig.tight_layout()
        plt.savefig("comPhase3647",dpi=720, bbox_inches='tight', pad_inches=0.2)

    def showRelativePhaseSingleIon(self, comFreInCm1, comFreInCm2, shift=0):
        fig= plt.figure(figsize=(8,3))
        bax = brokenaxes(xlims=((0, 50), (450,500), (900, 950)),ylims=((-2,1.5),))
        lab={'Ch0':'Mass1','Ch1':'Mass2','Ch2':'Mass16 x 10','Ch3':'Mass17 x 0.1','Ch4':'Mass18 x 0.1'}
        i=0
        for gas in self.spectraBottleB.keys():
            #if 'rebin' not in gas or 'filter' not in gas:
            if 'fil' not in gas or 'window' not in gas or 'rebin' in gas:
                continue
            elif 'Ch1' not in gas:
                continue
            label=lab[gas[14:]]
            if 'com' in gas:
                label = 'filter'+str(comFreInCm1)+'_'+label
                delay = self.delayB*comFreInCm1/comFreInCm2
                spec = self.spectraBottleB[gas]
            else:
                label = 'filter'+str(comFreInCm2)+'_'+label
                delay = self.delayB
                spec = self.spectraBottleB[gas]

            i=i+1                
            rebindelay = self.rebin_delay
            #else:
            #    delay = self.rebin_delay
            bax.plot(delay*10**15+shift, spec,label=label)
            bax.set_xlabel("Delay/fs")
            bax.set_ylabel('a.u.')
        bax.legend(ncol=2, loc=3)
        #plt.ylim([-1.5,1.5])
        #fig.tight_layout()
        plt.savefig("comPhase"+str(comFreInCm1)+'-'+str(comFreInCm2),dpi=720, bbox_inches='tight', pad_inches=0.2)

    def STFTS(self, gas, windowsize=0, ratio=1):
        '''
        windowsize is window size in fs 
        '''
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
            f, t, self.stftSB[gas] = sps.stft(self.spectraBottleB[gas], fs=1/(
                self.rebin_delay[2]-self.rebin_delay[1]), noverlap=windowsize-2, nperseg=windowsize, nfft=windowsize*3)
            t=t+self.windowSize*10**-15
            self.stftSBFre = f
            self.stftSBDelay = t
            vmax=abs(self.stftSB[gas][0:int(len(f)/10)]).max()*ratio
            vmin=abs(self.stftSB[gas]).min()*ratio
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            levels = np.arange(vmin,vmax,(vmax-vmin)/500)
            im=ax.contourf(t*10**15, f/10**12*33.35641,
                           np.abs(self.stftSB[gas]), levels=levels, cmap='jet', norm=norm)
            ax_histx.plot(t*10**15, np.sum(np.abs(self.stftSB[gas]),axis=0))
            ax_histx.set_xticklabels([])
            ax_histy.plot(np.sum(np.abs(self.stftSB[gas]),axis=1),f/10**12*33.35641)
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
            plt.savefig("sfts_180fs.png",dpi=720, bbox_inches='tight', pad_inches=0.2)
            plt.show()

    def cwtS(self, gas):
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
        scales = np.arange(20, 6000)/10
        #if 'rebin' not in gas or 'filter' not in gas:*
        if 'rebin' in gas:
            self.waveletsCoe[gas], freqs=pywt.cwt(self.spectraBottleB[gas],scales,'cmor10-0.5',sampling_period=self.rebin_delay[2]-self.rebin_delay[1])
            #print(self.waveletsCoe[gas].shape)
            t= self.rebin_delay
            vmax=abs(self.waveletsCoe[gas][:-5500]).max()
            vmin=abs(self.waveletsCoe[gas][:-5500]).min()
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            levels = np.arange(vmin,vmax,(vmax-vmin)/50)
            im=ax.contourf(t*10**15, freqs[:-5500]/10**12*33.35641,
                           np.abs(self.waveletsCoe[gas][:-5500]), levels=levels, cmap='jet', norm=norm)
            ax_histx.plot(t*10**15, np.sum(np.abs(self.waveletsCoe[gas]),axis=0))
            ax_histx.set_xticklabels([])
            ax_histy.plot(np.sum(np.abs(self.waveletsCoe[gas]),axis=1),freqs/10**12*33.35641)
            #plt.xlim([100, self.delayB[-1]])
            #ax_histy.set_ylim([200, 4500])
            ax.set_ylabel('Frequency [cm-1]')
            ax.set_xlabel('Time [fs]')
            #ax_histy.set_yticklabels([])
            ax_histy.yaxis.tick_right()
            #plt.clim(0,abs(self.stftSB[gas][:int(len(f)/15)]).max())
            #plt.colorbar(im, cax=a1, ticks=levels[::100])
            #plt.title(gas+d)
            fig.tight_layout()
            plt.savefig("wavelets_180fs.png",dpi=720, bbox_inches='tight', pad_inches=0.2)
            plt.show()

    def convolve_sft(self, gas, frequency=[3650], width=250):#https://stackoverflow.com/questions/24194238/how-can-i-find-the-fwhm-of-a-peak-in-a-noisy-data-set-in-python-numpy-scipy
        f_stfs = np.sum(np.abs(self.stftSB[gas]),axis=1)
        plt.figure(figsize=(5,5))
        for fre in frequency:
            if fre==1812:
                width = 900
                ratio = 1
            elif fre==4150:
                ratio = 1
                width = 250
            else:
                ratio = 1
                width = 350
            range_index=range((np.abs(self.stftSBFre/10**12*33.35641-(fre-width))).argmin(), (np.abs(self.stftSBFre/10**12*33.35641-(fre+width))).argmin(), 1) #https://stackoverflow.com/questions/8914491/finding-the-nearest-value-and-return-the-index-of-array-in-python
            f_stfs_i = f_stfs[range_index]
            FWHM=self.find_fwhm(f_stfs_i)
            sigma = FWHM/2*sqrt(2*log(2))

            convolve_window = np.reshape(sps.windows.gaussian(len(range_index), std=sigma), newshape=(len(range_index),1))
            convS = sps.convolve2d(np.abs(self.stftSB[gas])[range_index],convolve_window)
            plt.plot(f_stfs_i)
            plt.plot(convolve_window*30)
            #plt.plot(convS[:,0])
            plt.show()
            plt.plot(self.stftSBDelay*10**15,np.sum(convS,axis=0)*ratio, label='frequency:'+str(fre)+'x'+str(ratio))
        plt.legend()
        plt.xlabel('Delay/fs')
        plt.ylabel('a.u.')
        plt.savefig('sfts_lineout'+'180fs'+'.png',dpi=720, bbox_inches='tight', pad_inches=0.2)
        plt.show()

    def find_fwhm(self, peak):
        peakPOS = peak.argmax()
        difference = max(peak) - min(peak)
        HM = np.abs(difference / 2)
        nearestL = (np.abs(peak[:peakPOS] - HM)).argmin()
        nearestR = (np.abs(peak[peakPOS:] - HM)).argmin()
        FWHM=np.abs(nearestR-nearestL)
        return FWHM
    
    def phaseRetrive(self,gas,omega=[525.62772,805.21693,883.50191,1338.30036,1591.79458,1696.17455,2147.24515,2669.14501,3205.9563,3645.84333]):
        for gas in gas:
            interPhase=[]
            [f_window, Y_window, P_window] = self.fftSB['window_'+gas]
            f_window = f_window*33.35641
            amp = np.array([np.sum(np.where(np.logical_and(f_window>=ff-4, f_window<=ff+4),Y_window,0)) for ff in omega])
            #plt.plot(omega,amp)
            #plt.title("Fre_Amp")
            #plt.show()
            amp = np.reshape(amp,(1,len(amp)))
            osc = self.spectraBottleB_noPad['window_'+gas]
            hwindow = np.hanning(len(self.delayB_noPad))

            for i in range(len(self.delayB_noPad)):
                step=self.delayB[1]-self.delayB[0]
                oscA=np.sum(np.dot(amp,np.array([np.sin(2*np.pi*f*(self.delayB_noPad-self.delayB_noPad[0]+i*step)) for f in np.array(omega)*29979245800])),0)
                #sps.windows.flattop(len(window),sym=FALSE)
                #plt.plot(self.delayB_noPad,hwindow*(oscA/100))
                #plt.plot(self.delayB_noPad,osc)
                #plt.show()
                y,t=self.inter_padding(hwindow*(oscA/np.amax(oscA))*0.5+osc/np.amax(polynomial(osc, order=15, plot=False)), self.delayB_noPad, paddingSize = 20000)
                f_inter,Y,_ = self.FFT(t,y)
                #plt.plot(f_inter*33.35641,Y[:-1])
                #plt.xlim([300,4500])
                #plt.show()
                interPhase = interPhase+[Y]
            f_inter = f_inter*33.35641
            #plt.plot(f_inter)
            #plt.show()
            index0 = np.argmin(np.abs(f_inter-300))
            index1 = np.argmin(np.abs(f_inter-4500))
            self.phaseBottle = {}
            i=0
            if ifsavePhase:
                self.wks = op.new_sheet('w',lname=str('Phase')+str('_')+str(gas)+str('_')+str('%.1E' % Decimal(self.intensity))+str('_')+self.saveRef)
            for f in omega:
                i=i+1
                self.phaseBottle[str(f)] = np.sum(np.array(interPhase)[:,np.argmin(np.abs(f_inter-f+4)):np.argmin(np.abs(f_inter-f-4))],1)
            
            #print(index0,index1)
                if ifsavePhase:
                    self.wks.from_list(0, self.delayB_noPad*1e15, lname='Delay (fs)', axis='X')
                    self.wks.from_list(i, self.phaseBottle[str(f)]/np.amax(self.phaseBottle[str(f)]), lname=str(f), axis='Y')
                #instantPhase=np.angle(sps.hilbert(self.smooth2(self.phaseBottle[str(f)],5)))
                #aa=self.phaseBottle[str(f)]
                #plt.plot(aa)
                #plt.show()
            #return self.delayB_noPad, f_inter[index0:index1], np.array(interPhase)[:,index0:index1]
    def superlet(self,gas):
        for gas in gas:
            fs = 1/(self.delayB[1]-self.delayB[0])  # sampling frequency
            signal = self.spectraBottleB_noPad[gas]
            # frequencies of interest in Hz
            foi = np.linspace(int(300/33.35641)*1e12, int(4500/33.35641)*1e12, 300)
            scales = scale_from_period(1 / foi)

            spec = superlet(
                signal,
                samplerate=fs,
                scales=scales,
                order_max=100,
                order_min=1,
                c_1=0.8,
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

    def phaseRetrive2(self,omega=[3645.84333]):
        self.filterSpecMass1 = {}
        self.filterSpecMass17 = {}
        self.filterSpecMass18 = {}
        self.filterPhaseMass1 = {}
        self.filterPhaseMass17 = {}
        self.filterPhaseMass18 = {}
        for f in omega:
            f0=(f-0.1)/33.35641*1e12
            f1=(f+0.1)/33.35641*1e12
            self.filterSpecMass1[str(f)] = self.useFilter2(self.spectraBottleB_noPad["Ch0"],f0,f1)
            self.filterSpecMass17[str(f)] = self.useFilter2(self.spectraBottleB_noPad["Ch6"],f0,f1)
            self.filterSpecMass18[str(f)] = self.useFilter2(self.spectraBottleB_noPad["Ch8"],f0,f1)
            self.filterPhaseMass1[str(f)] = np.angle(sps.hilbert(self.filterSpecMass1[str(f)]))
            self.filterPhaseMass17[str(f)] = np.angle(sps.hilbert(self.filterSpecMass17[str(f)]))
            self.filterPhaseMass18[str(f)] = np.angle(sps.hilbert(self.filterSpecMass18[str(f)]))
            print(str(f)+": "+str(self.filterPhaseMass1[str(f)][1000]/np.pi)+str(", ")+str(self.filterPhaseMass17[str(f)][1000]/np.pi)+str(", ")+str(self.filterPhaseMass18[str(f)][1000]/np.pi))
            plt.plot(self.filterSpecMass1[str(f)])
            plt.plot(self.filterSpecMass17[str(f)])
            plt.plot(self.filterSpecMass18[str(f)])
            plt.show()

if __name__ == '__main__':
    
    f38 = [r'scan_tof_2022-06-04-22-37-09', r'scan_tof_2022-06-04-15-39-25']#r'scan_tof_2022-06-04-22-37-09', r'scan_tof_2022-06-04-15-39-25'
    f38r =[r'scan_tof_2022-06-03-22-44-54',r'scan_tof_2022-06-04-18-15-48',r'scan_tof_2022-06-01-20-53-44',r'scan_tof_2022-06-04-12-44-47']
    f38rr=[r'scan_tof_2022-06-04-15-39-25',r'scan_tof_2022-06-03-22-44-54',r'scan_tof_2022-06-04-18-15-48']
    f0 = [r'scan_tof_2022-07-21-20-43-45']#,r'scan_tof_2022-07-21-18-03-55',r'scan_tof_2022-07-22-09-58-27',r'scan_tof_2022-07-22-13-16-09']#150-150
    f1 = [r'scan_tof_2022-07-22-19-25-26']#,r'scan_tof_2022-07-22-16-06-00',r'scan_tof_2022-07-25-10-00-28']#180-180
    f2 = [r'scan_tof_2022-07-21-09-39-34']#,r'scan_tof_2022-07-21-12-38-03', r'scan_tof_2022-07-21-15-17-55']#200-200
    f3 = [r'scan_tof_2022-07-27-12-17-05 - Copy']

    f0 = FFT_ionS(filename=f0, scanTime=10, sampRate=300, molecule='H2', intensity=cal_intensity(150,13,5),ChN=11)
    f1 = FFT_ionS(filename=f1, scanTime=10, sampRate=300, molecule='H2', intensity=cal_intensity(180,11,5),ChN=11)
    f2 = FFT_ionS(filename=f2, scanTime=10, sampRate=300, molecule='H2', intensity=cal_intensity(200,10,5),ChN=11)
    f3 = FFT_ionS(filename=f3, scanTime=10, sampRate=300, molecule='H2', intensity=cal_intensity(230,9,5),ChN=11)

    #ff=[fp1,fp2,fp3,fp4]
    ff=[f2]
    for x in ff:
        x.pathFinder()
        print(x.filepath)
        x.read_splitB(isCorrection=False,sumN=1)
        x.checkData()
        #x.read_splitD()
        #x.mass_spectra()
        #plt.plot(np.arange(x.dataD.shape[0]),np.sum(x.dataD,1))
        #plt.plot(x.calculator.pixel2mass(np.arange(x.data.shape[0])),np.sum(x.data,1),label=None)
        #plt.show()
        #x.useFilter(100/33.35641*1e12, 6000/33.35641*1e12)

        
        x.window(windowSize=90, direction='left')
        x.rmvExp()
        #x.smooth(windowSize=9)
        #x.show_Spectra()
        #plt.show()
        
        #x.interpS(5)

        
        #x.rebinS(factor=5)
        #x.cwtS(gas='rebin_window_Ch8')
        #x.STFTS(gas='rebin_filter_Ch4', windowsize=int(9.14626*30/1200*2000), ratio=1)
        #x.convolve_sft(gas='rebin_filter_Ch4', frequency=[3650])
        
        
        #x.showRelativePhase()
        #x.showRelativePhaseSingleIon(comFreInCm1=4147,comFreInCm2=3647)
        #plt.show()
        
        #x.rmvExp()
        #
        x.padding(paddingSize=100000)
        x.FFTS()
        x.show_FFT()
        #x.phaseRetrive2()
        #x.superlet(['window_Ch8'])
        #plt.show()
        #x.phaseRetrive(['Ch10'],omega=[3647])
        #x.phaseRetrive(['Ch0','Ch6','Ch8'])
        #plt.plot(x.delayB_noPad, x.phaseBottle['3647'])
        #plt.show()
        #x.show_FFTD()
        #plt.show()
        
        #x.show_FFT_comp()
        #op.save(pl.PureWindowsPath(x.exportPath,r'20220602H2O.opju'))
        #            if ifsavePhase:
        #        self.wk = op.new_book('w',lname=str('Phase')+str(gas)+str('%.1E' % Decimal(self.intensity)))
        #    for f in omega:
        #        self.phaseBottle[str(f)] = np.sum(np.array(interPhase)[:,np.argmin(np.abs(f_inter-f+20)):np.argmin(np.abs(f_inter-f-20))],1)
        #    #print(index0,index1)
        #        if ifsavePhase:
        #            self.wks = self.wk.add_sheet(name=str(f))
        #            self.wks.from_list(0, self.delayB_noPad*1e15, lname='Delay (fs)', axis='X')
        #            self.wks.from_list(1, self.phaseBottle[str(f)]/np.amax(self.phaseBottle[str(f)]), lname=str(f), axis='Y')
        
    