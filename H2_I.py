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
from decimal import Decimal
from cal_intensity import cal_intensity
from calculate_k_b import Calibration_mass

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
    def __init__(self, filename, scanTime, sampRate, massRange=[5, 50], molecule='H2', intensity = 0, dcRange=5, ChN=11, cal_mass=[2,16], cal_pixel=[3082,10216]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[494,583]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[197,279]):
        '''

        '''
        self.filename = filename
        self.saveRef = str(filename)
        self.filepath = []
        self.scanTime = scanTime
        self.trueScanTime = ceil(self.scanTime/1.6)*1.6
        self.scanLengthB = 3360
        #self.delayB = np.arange(self.scanLengthB)/self.scanLengthB*100*2*2*3.33564*10**-15/4161.16632*4155.0587
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
        self.rootPath= pl.PureWindowsPath(r'D:\SF_FFT')
        self.savePath = os.path.join(r'D:\SF_FFT', r'processedData', r'inter', self.molecule, str(self.intensity)+'_'+self.molecule)

        
    def pathFinder(self):
        for fn in self.filename:
            interpath = pl.PureWindowsPath(fn[9:16], fn[9:19])
            self.filepath = self.filepath + [pl.PureWindowsPath(self.rootPath,r'OriginalData',interpath, fn + r'.hdf5')]
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

    def read_splitB(self, isCorrection = True, sumN=1):
        ChN=self.ChN #number of channels measured
        for fp in self.filepath:
            self.delayB = np.arange(self.scanLengthB)/self.scanLengthB*100*2*2*3.33564*10**-15/4161.16632*4146.12
            if isCorrection:
                with h5py.File(fp, 'r+') as f:
                    print(f.keys())
                    data = np.array(f['dataB'], dtype=float)
                    print(data.shape)
                    m = int(data.shape[0]/self.ChN/sumN)
                    for i in [0,2,4,6,8,10]:
                        self.interSpectraBottleB['Ch'+str(i)] = data[i::ChN][-int(m*sumN):].reshape(m, sumN, data.shape[1]).sum(axis=1)
                    self.interinterSpectraBottleB = {}
                    num = int(m/4)
                    for gas in self.interSpectraBottleB.keys():
                        self.interinterSpectraBottleB[gas] = {}
                        self.spectraBottleB[gas] = np.zeros((num, 13000))
                    for i in range(num):
                        for gas in self.interSpectraBottleB.keys():
                            _inter = self.interSpectraBottleB[gas][i*4:(i+1)*4]
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
        '''
        calibrate by time zero
        '''
        plt.clf()
        _iS = np.array(np.zeros(self.interSpectraBottleB[gas].shape))
        #_iS2 = np.array(np.zeros(self.interSpectraBottleB[gas].shape))
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            #_iS2[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], (4161.07887-20)/33.35641*1e12, (4161.07887+20)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            _iS[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], 100/33.35641*1e12, 4000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
        #iS2 = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS2)
        _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=20000)
        indexMax = []
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange)
            indexMax = indexMax + [np.argmax(_inter)]
            #plt.plot(_delayRange,_inter)
            
        print(indexMax)
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        #plt.show()
        _ref = sum(indexMax[int(self.interSpectraBottleB[gas].shape[0]/2)-5:int(self.interSpectraBottleB[gas].shape[0]/2)+5])/10
        _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/20000
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange+_shift[i])
            #delayRangeA = np.linspace(start=200,stop=1000, num=50000)#Used to compare the zreo delay when calibrating using stretch mode frequency to calibrate
            #_inter2 = iS2.ev(i,delayRangeA+_shift[i])#
            #plt.plot(delayRangeA,_inter2)
            #plt.plot(_delayRange,_inter)
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        #plt.show()
        return _shift

    def calDrift2(self, _cRange, gas='Ch0'):
        '''
        calibrate by drift of the strech mode oscillation
        '''
        plt.clf()
        _iS = np.array(np.zeros(self.interSpectraBottleB[gas].shape))
        _iS2 = np.array(np.zeros(self.interSpectraBottleB[gas].shape))
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _iS[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], (10)/33.35641*1e12, (5000)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            #_iS[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], (4161.07887-200)/33.35641*1e12, (4161.07887+200)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            _iS2[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], 100/33.35641*1e12, 4000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
        iS2 = sp.interpolate.RectBivariateSpline(range(_iS2.shape[0]), self.delayB*1e15, _iS2)
        _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=2000)
        indexMax = []
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange)
            plt.plot(_delayRange,_inter)
            #indexMax = indexMax + [sps.argrelextrema(_inter, np.greater)[0][-1]]
            ref1=np.argmax(np.abs(_inter))
            aveRange=100
            ref2=np.dot(_inter[ref1-aveRange:ref1+aveRange],np.arange(-aveRange,aveRange,1))/np.sum(_inter[ref1-aveRange:ref1+aveRange])
            indexMax = indexMax + [int(ref1+ref2)]
            
        print(indexMax)
        plt.xlabel('Delay (fs)')
        plt.ylabel('a.u.')
        #plt.show()
        _ref = sum(indexMax[int(self.interSpectraBottleB[gas].shape[0]/2)-5:int(self.interSpectraBottleB[gas].shape[0]/2)+5])/10
        _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/2000
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange+_shift[i])
            delayRangeA = np.linspace(start=400,stop=500, num=2000)#Used to compare the zreo delay when calibrating using stretch mode frequency to calibrate
            _inter2 = iS2.ev(i,delayRangeA+_shift[i])#
            #plt.plot(_delayRange,_inter)
            plt.plot(delayRangeA,_inter2)
        plt.xlabel('Delay (fs)')
        plt.ylabel('a.u.')
        #plt.show()
        return _shift
    
    def calDrift3(self, _cRange, gas='Ch0'):
        '''
        calibrate by drift of the strech mode oscillation
        '''
        op.set_show(show=True)
        plt.clf()
        _iS = np.array(np.zeros(self.interSpectraBottleB[gas].shape))
        _iS2 = np.array(np.zeros(self.interSpectraBottleB[gas].shape))
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _iS[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], 500/33.35641*1e12, 5000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            #_iS[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], (4161.07887-200)/33.35641*1e12, (4161.07887+200)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            #_iS2[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], 100/33.35641*1e12, 4000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
        #iS2 = sp.interpolate.RectBivariateSpline(range(_iS2.shape[0]), self.delayB*1e15, _iS2)
        _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=2000)
        indexMax = []
        self.wks = op.new_sheet('w',lname='original'+str('_')+str(self.printNum))
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange)
            plt.plot(_delayRange,_inter)
            #indexMax = indexMax + [sps.argrelextrema(_inter, np.greater)[0][-1]]
            ref1=np.argmax(np.abs(_inter))
            aveRange=50
            ref2=np.dot(_inter[ref1-aveRange:ref1+aveRange],np.arange(-aveRange,aveRange,1))/np.sum(_inter[ref1-aveRange:ref1+aveRange])
            indexMax = indexMax + [int(ref1+ref2)]
            
            self.wks.from_list(0, _delayRange, 'X')
            self.wks.from_list(i+1, _inter, lname=str(i+1), axis='Y')
        print(indexMax)
        plt.xlabel('Delay (fs)')
        plt.ylabel('a.u.')
        plt.show()
        _ref = sum(indexMax[int(self.interSpectraBottleB[gas].shape[0]/2)-5:int(self.interSpectraBottleB[gas].shape[0]/2)+5])/10
        _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/2000
        self.wks = op.new_sheet('w',lname='calibrated'+str('_')+str(self.printNum))
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange+_shift[i])
            #delayRangeA = np.linspace(start=50,stop=100, num=2000)#Used to compare the zreo delay when calibrating using stretch mode frequency to calibrate
            #_inter2 = iS2.ev(i,delayRangeA+_shift[i])#
            plt.plot(_delayRange,_inter)
            #plt.plot(delayRangeA,_inter2)
            
            self.wks.from_list(0, _delayRange, 'X')
            self.wks.from_list(i+1, _inter, lname=str(i+1), axis='Y')
        plt.xlabel('Delay (fs)')
        plt.ylabel('a.u.')
        plt.show()
        return _shift

    #def delayCorrection(self, _cRange = [400,404]):
    def delayCorrection(self, _cRange = [20,100]):
        print('self.num is'+str(self.num)+'!/n')
        for k in range(self.num):
            for gas in self.interSpectraBottleB.keys():
                self.interSpectraBottleB[gas] = self.interinterSpectraBottleB[gas][str(k)]
            xxx = 0
            while True:
                try:
                    print(_cRange)
                    #_shift = self.calDrift(_cRange = _cRange)
                    _shift = self.calDrift2(_cRange = _cRange)
                    break
                except IndexError:
                    print('Wrong Region! Try again!')
                    if xxx>10:
                        _cRange = [np.array(_cRange)[0], np.array(_cRange)[1]+0.5]
                    _cRange = np.array(_cRange)+0.5
                    xxx = xxx+1
            #############################################################
            _delayRange = np.linspace(start=0,stop=1300, num=13000)# New delay!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ###############################################################
            for gas in self.interSpectraBottleB.keys():
                iS = sp.interpolate.RectBivariateSpline(range(self.interSpectraBottleB[gas].shape[0]), self.delayB*1e15, self.interSpectraBottleB[gas])
                inter = 0
                for i in range(self.interSpectraBottleB[gas].shape[0]):
                    inter = inter + iS.ev(i, _delayRange+_shift[i])
                if gas == 'Ch2':
                    print(sum(inter))
                self.spectraBottleB[gas][k] = inter

        self.delayB = _delayRange * 1e-15
        self.interSpectraBottleB = self.spectraBottleB.copy()

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

if __name__ == '__main__':
    
    f0 = [r'scan_tof_2022-07-21-20-43-45',r'scan_tof_2022-07-21-18-03-55',r'scan_tof_2022-07-22-09-58-27',r'scan_tof_2022-07-22-13-16-09']#150-150
    f1 = [r'scan_tof_2022-07-22-16-06-00',r'scan_tof_2022-07-22-19-25-26',r'scan_tof_2022-07-25-10-00-28']#180-180
    f2 = [r'scan_tof_2022-07-21-09-39-34',r'scan_tof_2022-07-21-12-38-03', r'scan_tof_2022-07-21-15-17-55']#200-200

    f0=FFT_ionS(filename=f0, scanTime=10, sampRate=300, molecule='H2', intensity=cal_intensity(150,13,5),ChN=11)
    f1=FFT_ionS(filename=f1, scanTime=10, sampRate=300, molecule='H2', intensity=cal_intensity(180,11,5),ChN=11)
    f2=FFT_ionS(filename=f2, scanTime=10, sampRate=300, molecule='H2', intensity=cal_intensity(200,10,5),ChN=11)

    #ff=[fp1,fp2,fp3,fp4]
    ff=[f0,f1,f2]
    for x in ff:
        x.pathFinder()
        print(x.filepath)
        
        x.read_splitB(isCorrection=True,sumN=3)
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
        
    