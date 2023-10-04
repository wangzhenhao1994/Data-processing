########################################################################
#This script is used to read in original data, then calibrate the delay drift for the first time, 'H2O_II.py' will calibrate again.
#ignore the functions not called.
#The processed files are saved in the folder 'processedData\inter\H2O', which will then used by script 'H2O_II.py'
#
#########################################################################
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
    #def __init__(self, filename, massRange=[5, 90], cal_mass=[45,61], cal_pixel=[415,1141]):
    def __init__(self, filename, scanTime, sampRate, massRange=[5, 50], molecule='H2O', intensity = 0, dcRange=5, ChN=11, cal_mass=[2,16], cal_pixel=[3082,10216]):
    #def __init__(self, filename, massRange=[5, 50], cal_mass=[17,18], cal_pixel=[494,583]):
    #def __init__(self, filename, massRange=[5, 50], cal_mass=[17,18], cal_pixel=[197,279]):
        '''
        scanTime: the scanTime set in the measurement program. It is not true scan time because the scan time can only be integral times of 1.6 second.
        sampRate: sample rate of the NIDaq.
        intensity: int or [a,b] when pump and probe intensity are different
        massRange=[5, 50], cal_mass=[2,16], cal_pixel=[3082,10216]: only used when reading data of digitizer.
        '''
        self.filename = filename
        self.saveRef = str(filename)
        self.filepath = []
        self.scanTime = scanTime
        self.trueScanTime = ceil(self.scanTime/1.6)*1.6
        self.scanLengthB = int(self.trueScanTime*sampRate)
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

        #setting for data from digitizer
        #self.calculator = Calibration_mass(cal_mass, cal_pixel)
        #self.channelSize = 12032#24000#1536
        #self.scanLengthD = 3320#1200
        #self.peakRange = [-100, 100]  # range of the peak
        #self.delayD = np.arange(self.scanLengthD)/self.scanLengthD*100*2*2*3.33564*10**-15
        #self.gasBottle = {
        #    "O+": self.calculator.cal_pixel(16)+self.peakRange,
        #    "H2+": self.calculator.cal_pixel(2)+self.peakRange,
        #    "H+": self.calculator.cal_pixel(0.99)+self.peakRange,
        #    "O+H2": #[0,self.channelSize]
        #    np.append(self.calculator.cal_pixel(16) + self.peakRange, self.calculator.cal_pixel(2) + self.peakRange),
        #    "O-H2": #[0,self.channelSize]
        #    np.append(np.append(self.calculator.cal_pixel(16) + self.peakRange, self.calculator.cal_pixel(2) + self.peakRange), 0),
        #}
        #self.spectraBottleD = {}
        #self.fftSD = {}
        #self.stftSD = {}
        #self.dataD = 0

        self.rootPath= pl.PureWindowsPath(r'D:\SF_FFT')
        self.savePath = os.path.join(self.rootPath, r'processedData', r'inter', self.molecule, str(self.intensity)+'_'+self.molecule)

    def pathFinder(self):
        for fn in self.filename:
            interpath = pl.PureWindowsPath(fn[9:16], fn[9:19])
            self.filepath = self.filepath + [pl.PureWindowsPath(self.rootPath, r'OriginalData', interpath, fn + r'.hdf5')]

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
        '''
        For data measured using boxcar.
        '''
        ChN=self.ChN #number of channels measured
        for fp in self.filepath:
            self.delayB = np.arange(self.scanLengthB)/self.scanLengthB*100*2*2*3.33564/4161.16632*4146.12#Change the delay so that the frequency of H2 foundamental vibration is 4146.12.
            self.stepSize = np.mean(self.delayB[1:]-self.delayB[:-1])
            if isCorrection:
                with h5py.File(fp, 'r+') as f:
                    print(f.keys())
                    data = np.array(f['dataB'], dtype=float)
                    print(data.shape)
                    #sumN: every sumN scans are add together before delay calibration.
                    #self.ChN: In total the data of 12 analog input channels are saved. Channels 0,2,4,6,8 are data of Mass1, Mass2, Mass16 and Mass17. Channel 10 is the reference data for Mass2.
                    #Channels 1,3,5,7,9,11 are grounded analog input channels inbetween every two channels for signal from boxcars, used to make sure the sample rate is low enough.
                    m = int(data.shape[0]/self.ChN/sumN)
                    for i in [0,2,4,6,8,10]:
                        self.interSpectraBottleB['Ch'+str(i)] = data[i::ChN][-int(m*sumN):].reshape(m, sumN, data.shape[1]).sum(axis=1)
                    self.interinterSpectraBottleB = {}
                    num = int(m/10)#every 59*sumN scans are delay-calibrated.
                    for gas in self.interSpectraBottleB.keys():
                        self.interinterSpectraBottleB[gas] = {}
                        self.spectraBottleB[gas] = np.zeros((num, 13000))
                    for i in range(num):
                        for gas in self.interSpectraBottleB.keys():
                            _inter = self.interSpectraBottleB[gas][i*10:(i+1)*10]
                            self.interinterSpectraBottleB[gas][str(i)] = _inter 
            self.num = num
            self.delayCorrection(file=os.path.split(fp)[1])
            if not os.path.exists(self.savePath):
                os.mkdir(self.savePath)
            print(self.interSpectraBottleB['Ch0'].shape)
            save_obj(self.interSpectraBottleB, pl.PureWindowsPath(self.savePath, os.path.split(fp)[1].replace(r'.hdf5','')+r'.pkl'))

    def read_splitD(self, overNight = True, firstTry = False, sumNo = 89, usefulRange = [0,2], cu = False):
        '''
        For data measured using digitizer.
        '''
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
        '''
        For data measured using digitizer.
        '''
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

    def calDrift(self, file, _cRange, gas='Ch0'):
        '''
        calibrate by drift of the strech mode oscillation
        '''
        plt.clf()
        calNum = 20000
        standard = 0
        standard = self.interSpectraBottleB['Ch8']-self.interSpectraBottleB['Ch0']#-self.interSpectraBottleB['Ch6']
        _iS = np.array(np.zeros(self.interSpectraBottleB[gas].shape))
        #_iS2 = np.array(np.zeros(self.interSpectraBottleB[gas].shape))
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            #_iS[i]=self.butter_bandpass_filter(standard[i], (3655.52723-40)/33.35641*1e12, (3655.52723+40)/33.35641*1e12, 1/self.stepSize*1e15)
            #_iS[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], (3655.52723-30)/33.35641*1e12, (3655.52723+30)/33.35641*1e12, 1/self.stepSize*1e15)
            _iS[i]=self.butter_bandpass_filter(standard[i], 100/33.35641*1e12, 1500/33.35641*1e12, 1/self.stepSize*1e15)
            #_iS2[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], 100/33.35641*1e12, 4000/33.35641*1e12, 1/self.stepSize)
        iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB, _iS)
        #iS2 = sp.interpolate.RectBivariateSpline(range(_iS2.shape[0]), self.delayB, _iS2)
        _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=calNum)
        indexMax = []
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange)
            #plt.plot(_delayRange,_inter)
            #indexMax = indexMax + [sps.argrelextrema(_inter, np.greater)[0][-1]]
            indexMax = indexMax + [np.argmax(np.abs(_inter))]
            
        print(indexMax)
        #plt.title('Before calibration')
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        #plt.show()
        _ref = np.mean(indexMax)
        _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/calNum
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange+_shift[i])
            #delayRangeA = np.linspace(start=50,stop=100, num=2000)#Used to compare the zreo delay when calibrating using stretch mode frequency to calibrate
            #_inter2 = iS2.ev(i,delayRangeA+_shift[i])#
            plt.plot(_delayRange,_inter)
            #plt.plot(delayRangeA,_inter2)
        plt.title('After calibration')
        plt.xlabel('Delay (fs)')
        plt.ylabel('a.u.')
        plt.savefig(pl.PureWindowsPath(self.savePath, os.path.split(file)[1].replace(r'.hdf5','')+r'.png'))
        #plt.show()
        return _shift

    #def delayCorrection(self, file, _cRange = [300,307]):#using frequency of sym.Strech as standard
    def delayCorrection(self, file, _cRange = [1,150]):#using auto correlation as standard
        for k in range(self.num):
            for gas in self.interSpectraBottleB.keys():
                self.interSpectraBottleB[gas] = self.interinterSpectraBottleB[gas][str(k)]
            xxx = 0
            while True:
                try:
                    print(_cRange)
                    #_shift = self.calDrift(_cRange = _cRange)
                    _shift = self.calDrift(file=file, _cRange = _cRange)
                    break
                except IndexError:
                    print('Wrong Region! Try again!')
                    if xxx>10:
                        _cRange = [np.array(_cRange)[0], np.array(_cRange)[1]+0.5]
                    _cRange = np.array(_cRange)+0.5
                    xxx = xxx+1
            #############################################################
            _delayRange,self.stepSize = np.linspace(start=0,stop=1300, num=13000,endpoint=False,retstep=True)# New delay!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ###############################################################
            for gas in self.interSpectraBottleB.keys():
                iS = sp.interpolate.RectBivariateSpline(range(self.interSpectraBottleB[gas].shape[0]), self.delayB, self.interSpectraBottleB[gas])
                inter = 0
                for i in range(self.interSpectraBottleB[gas].shape[0]):
                    inter = inter + iS.ev(i, _delayRange+_shift[i])
                if gas == 'Ch8':
                    print(sum(inter))
                self.spectraBottleB[gas][k] = inter
        self.delayB = _delayRange
        self.interSpectraBottleB = self.spectraBottleB.copy()

    def findZeroDelay(self):
        zeroIndex = {}
        zeroIndex2=[]
        for gas in list(self.interSpectraBottleB.keys()):
            if gas not in ['Ch0', 'Ch2', 'Ch4', 'Ch6', 'Ch8', 'Ch10']:
                continue
            _Spec = self.interSpectraBottleB[gas]
            _Spec=self.butter_bandpass_filter(_Spec, 1/33.35641*1e12, 150/33.35641*1e12, 1/self.stepSize)
            #plt.plot(_Spec)
            #plt.show()
            if gas in ['Ch0', 'Ch2', 'Ch4', 'Ch6']:
                zeroIndex[gas] =(-1*_Spec[:500]).argmax()
                zeroIndex2 = zeroIndex2 + [(-1*_Spec[:500]).argmax()]
            else:
                zeroIndex[gas] =_Spec[:500].argmax()
                zeroIndex2 = zeroIndex2 + [_Spec[:500].argmax()]
        return zeroIndex, zeroIndex2

    def show_Spectra(self, ifsaveT = False):
        gs = gridspec.GridSpec(2, 3)
        #gs = gridspec.GridSpec(1, 1)
        fig = plt.figure(figsize=(20,8))
        #ax = fig.add_subplot(111)
        lab=['Mass1','Mass2','Mass16','Mass17','Mass18','ref2']
        i=0
        if ifsaveT:
            op.set_show(show=True)
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
            ax.plot(delay*10**15 ,
                     self.spectraBottleB[gas].sum(axis=0), label=gas)#/np.amax(self.spectraBottleB[gas])
            #ax.plot(delay*10**15 ,
            #         self.spectraBottleB["window_"+gas], label=gas)#/np.amax(self.spectraBottleB[gas])
            #ax.plot(delay*10**15 ,
            #         self.spectraBottleB['window_filter_'+gas], label=label)
            #ax.plot(delay*10**15 ,
            #         self.spectraBottleB['window_'+gas], label='window_'+gas)
            #ax.plot(rebindelay*10**15 ,
            #         self.spectraBottleB['rebin_window_'+gas], label='rebin_window_'+gas)
            ax.set_xlabel("Delay/fs")
            ax.set_ylabel('a.u.')
            #plt.xlim([300,400])
            ax.legend()
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('Time')+label+self.molecule+self.stage+str('%.1E' % Decimal(self.intensity))+str('.dat')), x.spectraBottleB[gas])
            #np.savetxt(pl.PureWindowsPath(self.exportPath, str('delay')+str('%.1E' % Decimal(self.intensity))+str('.dat')), delay*10**15)
            if ifsaveT:
                self.wks.from_list(0, delay, 'X')
                self.wks.from_list(i, self.spectraBottleB[gas].sum(axis=0), lname=gas, axis='Y')
        #plt.legend()
        fig.tight_layout()
        #plt.savefig("spectra_31"+d,dpi=720, bbox_inches='tight', pad_inches=0.2)
        plt.show()

if __name__ == '__main__':
    
    #f0 = [r'scan_tof_2022-05-24-17-38-18']
    #f1=[r'scan_tof_2022-05-25-09-47-24']
    #f2=[r'scan_tof_2022-05-25-12-31-23']
    #f3=[r'scan_tof_2022-05-25-14-08-02']
    #f4=[r'scan_tof_2022-05-25-15-50-01']
    #f5=[r'scan_tof_2022-05-25-17-27-19']
    #f6=[r'scan_tof_2022-05-25-18-56-51']
    #f7=[r'scan_tof_2022-05-26-12-44-27']
    #f8=[r'scan_tof_2022-05-26-14-13-58']
    #f9=[r'scan_tof_2022-05-26-15-46-15']
    #f10=[r'scan_tof_2022-05-26-17-00-01',r'scan_tof_2022-05-26-22-52-13']
    f11=[r'scan_tof_2022-05-26-18-50-20']
    #f12=[]
    #f13=[r'scan_tof_2022-05-27-10-23-47']
    #f14=[r'scan_tof_2022-05-27-11-54-12',r'scan_tof_2022-05-27-19-53-05']
    f15=[r'scan_tof_2022-05-27-13-12-13']
    #f16=[r'scan_tof_2022-05-27-15-11-46']
    #f17=[r'scan_tof_2022-05-27-16-32-58']
    f18=[r'scan_tof_2022-05-27-18-03-09']
#
    ##since 20220530
    #
    #f19 = [r'scan_tof_2022-05-30-10-43-15']
    #f20 = [r'scan_tof_2022-05-30-12-23-32']
    #f21 = [r'scan_tof_2022-05-30-13-50-56']
    #f22 = [r'scan_tof_2022-05-30-15-27-08']
    #f23 = [r'scan_tof_2022-05-30-16-57-23']
    #f24 = [r'scan_tof_2022-05-30-18-18-33']
    #f25 = [r'scan_tof_2022-05-30-19-44-25']
    #f26 = [r'scan_tof_2022-05-31-09-23-21']
    #f27 = [r'scan_tof_2022-05-31-11-12-23']
    #f28 = [r'scan_tof_2022-05-31-12-26-11']    
    #f29 = [r'scan_tof_2022-05-31-14-16-24']
    #f30 = [r'scan_tof_2022-05-31-15-42-05']
    #f31 = [r'scan_tof_2022-05-31-17-51-43']
    #f32 = [r'scan_tof_2022-05-31-19-41-25']
    ##1st June, today the fiber output is a little bit less, and the data also shows lower intensity though try to set as same parameter as last round.
    f33 = [r'scan_tof_2022-06-01-10-40-05']
    #f34 = [r'scan_tof_2022-06-01-12-10-05']
    #f35 = [r'scan_tof_2022-06-01-13-39-07']
    #f36 = [r'scan_tof_2022-06-01-15-20-40']
    #f37 = [r'scan_tof_2022-06-01-16-54-40']
    #r'scan_tof_2022-06-04-22-37-09',r'scan_tof_2022-06-01-18-10-15',r'scan_tof_2022-06-04-12-44-47',r'scan_tof_2022-06-01-15-20-40',r'scan_tof_2022-06-01-20-53-44',r'scan_tof_2022-06-02-17-31-08',
    f38 = [r'scan_tof_2022-06-04-22-37-09', r'scan_tof_2022-06-04-15-39-25']#r'scan_tof_2022-06-04-22-37-09', r'scan_tof_2022-06-04-15-39-25'
    f38r =[r'scan_tof_2022-06-03-22-44-54',r'scan_tof_2022-06-04-18-15-48',r'scan_tof_2022-06-01-20-53-44',r'scan_tof_2022-06-04-12-44-47']
    #f38rr=[r'scan_tof_2022-06-04-15-39-25',r'scan_tof_2022-06-03-22-44-54',r'scan_tof_2022-06-04-18-15-48']
    #f38rrr = [r'scan_tof_2022-06-01-18-10-15',r'scan_tof_2022-06-02-18-58-56',r'scan_tof_2022-06-03-09-31-11',r'scan_tof_2022-06-03-12-20-18',r'scan_tof_2022-06-03-18-10-43',r'scan_tof_2022-06-01-15-20-40',r'scan_tof_2022-06-02-17-31-08']
    #r'scan_tof_2022-06-04-22-37-09',r'scan_tof_2022-06-01-18-10-15',r'scan_tof_2022-06-02-18-58-56',r'scan_tof_2022-06-03-09-31-11',,r'scan_tof_2022-06-03-12-20-18',r'scan_tof_2022-06-03-18-10-43',r'scan_tof_2022-06-01-15-20-40',r'scan_tof_2022-06-02-17-31-08'
    #f39 = [r'scan_tof_2022-06-02-10-12-43']
    #f43 = [r'scan_tof_2022-06-02-17-31-08']
    #f44 = [r'scan_tof_2022-06-02-18-58-56']
    #f45 = [r'scan_tof_2022-06-03-09-31-11',r'scan_tof_2022-06-03-10-56-37']
    #f46 = [r'scan_tof_2022-06-03-16-24-18']

    #f47 = []
    f48 = [r'scan_tof_2022-07-01-07-53-13',r'scan_tof_2022-06-21-18-10-59',r'scan_tof_2022-06-22-19-24-23',r'scan_tof_2022-06-23-10-09-03',r'scan_tof_2022-06-23-11-48-23',r'scan_tof_2022-06-27-19-55-38',r'scan_tof_2022-06-28-15-47-53',r'scan_tof_2022-06-28-19-18-27',r'scan_tof_2022-06-29-12-09-27',r'scan_tof_2022-06-29-15-28-34',r'scan_tof_2022-07-01-07-53-13'] #200-200
#
#
#
#

    #r'scan_tof_2022-06-06-18-27-03',r'scan_tof_2022-06-05-17-03-06', 200-200 #200-200
    #r'scan_tof_2022-06-07-09-56-51' 200-200 delay drift ...try to relocate the zero delay later
    #scan_tof_2022-06-07-11-29-36 200-200 delay little drift ...try to relocate the zero delay later
    #scan_tof_2022-06-07-13-33-30 250-230
    #f49 = [r'scan_tof_2022-06-07-16-29-08']
    f49 = [r'scan_tof_2022-06-07-13-33-30',r'scan_tof_2022-06-07-16-29-08',r'scan_tof_2022-06-07-17-41-15']#250-230
    #ftest = [r'scan_tof_2022-06-22-19-24-23']
    #250-230 good scan_tof_2022-06-07-13-33-30 scan_tof_2022-06-07-16-29-08
    ##250-230 bad scan_tof_2022-06-07-17-41-15

    #phaseRetrive
    #fp1=[r'scan_tof_2022-06-04-22-37-09']
    #fp2=[r'scan_tof_2022-06-04-15-39-25',r'scan_tof_2022-06-03-22-44-54',]
    #fp3=[r'scan_tof_2022-06-04-18-15-48']
    #fp4 = fp1+fp2+fp3


    ############################20220621
    #f50 = [r'scan_tof_2022-06-23-10-09-03',r'scan_tof_2022-06-23-11-48-23']
    #f51 = [r'scan_tof_2022-06-23-13-11-05'] #190-200
    #f52 = [r'scan_tof_2022-06-23-15-25-48'] #180-200
    f53 = [r'scan_tof_2022-06-05-13-10-09'
            ,r'scan_tof_2022-06-05-15-07-15'
            ,r'scan_tof_2022-06-23-17-22-11'
            ,r'scan_tof_2022-06-23-20-24-23'
            ,r'scan_tof_2022-06-24-11-11-20'
            ,r'scan_tof_2022-06-24-15-11-10'
            ,r'scan_tof_2022-06-24-20-31-22'
            ,r'scan_tof_2022-06-27-11-08-21'
            ,r'scan_tof_2022-06-27-13-40-24'
            ,r'scan_tof_2022-06-27-17-09-33'] #180-180
    #f54 = [r'scan_tof_2022-06-27-19-55-38']
    #f54 = [r'scan_tof_2022-06-28-12-45-25']
    #f54 = [r'scan_tof_2022-06-28-15-47-53']
    ##f54 = [r'scan_tof_2022-06-28-19-18-27']
    ##f54 = [r'scan_tof_2022-06-29-12-09-27']
    #f55 = [r'scan_tof_2022-06-29-18-16-59']
    #f56 = [r'scan_tof_2022-07-03-12-09-43']
    
    f0=[r'scan_tof_2022-08-25-11-06-12',r'scan_tof_2022-08-25-13-35-09',r'scan_tof_2022-08-25-16-25-51',r'scan_tof_2022-08-25-19-06-27',r'scan_tof_2022-08-25-21-44-21',r'scan_tof_2022-08-26-00-33-41']
    f1=[r'scan_tof_2022-08-26-08-42-21',r'scan_tof_2022-08-26-11-24-19',r'scan_tof_2022-08-26-14-10-20',r'scan_tof_2022-08-26-17-47-46',r'scan_tof_2022-08-26-21-02-27']
    f2=[r'scan_tof_2022-08-27-09-37-01',r'scan_tof_2022-08-27-15-02-52',r'scan_tof_2022-08-27-18-07-19',r'scan_tof_2022-08-27-21-01-19']
    ####################################################################################################################202209
    f00 = [r'scan_tof_2022-09-09-14-16-23',r'scan_tof_2022-09-11-11-29-48',r'scan_tof_2022-09-11-15-42-18']
    
    
    
   


    #fp1=FFT_ionS(filename=fp1, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    #fp2=FFT_ionS(filename=fp2, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    #fp3=FFT_ionS(filename=fp3, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    #fp4=FFT_ionS(filename=fp4, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
#
    #fall = f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18+f19+f20+f21+f22+f23+f24+f25+f26+f27+f28+f29+f30+f31+f32+f33+f34+f35+f36+f37+f38+f39+f43+f44+f45+f46+ftest
    #
    ##probe 230 mW
    #f0 = FFT_ionS(filename=f0, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)
    #f1 = FFT_ionS(filename=f1, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)
    #f2 = FFT_ionS(filename=f2, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(60,8.2,5),ChN=11)
    #f3 = FFT_ionS(filename=f3, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11)
    #f4 = FFT_ionS(filename=f4, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)
    #f5 = FFT_ionS(filename=f5, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)#repeat
#
    #f19 = FFT_ionS(filename=f19, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)
    #f20 = FFT_ionS(filename=f20, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.5,5),ChN=11)#maybe repeat again
    #f21 = FFT_ionS(filename=f21, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.5,5),ChN=11)
    #f22 = FFT_ionS(filename=f22, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.5,5),ChN=11)
#
    ##probe 100 mW
    #f6 = FFT_ionS(filename=f6, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)
    #f16 = FFT_ionS(filename=f16, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)
    #f17 = FFT_ionS(filename=f17, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)#repeat
    f18 = FFT_ionS(filename=f18, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11)
#
#
    ##probe 150 mW
    #f7 = FFT_ionS(filename=f7, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)
    #f13 = FFT_ionS(filename=f13, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)
    #f14 = FFT_ionS(filename=f14, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.5,5),ChN=11)#repeat
    f15 = FFT_ionS(filename=f15, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11)
    #
    #f34 = FFT_ionS(filename=f34, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)
    #f35 = FFT_ionS(filename=f35, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9.8,5),ChN=11)#repeat 20220601 actual intensity seems lower due to longer pulse
    #f36 = FFT_ionS(filename=f36, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)#good repeat 20220601
    #f37 = FFT_ionS(filename=f37, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,11.6,5),ChN=11)
    ##f38 = FFT_ionS(filename=f38, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    #f45 = FFT_ionS(filename=f45, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11) #actua power 155-155
    #
    ##probe 180mW
    #f8 = FFT_ionS(filename=f8, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11)#repeat
    #f9 = FFT_ionS(filename=f9, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)#interesting, multi fre in H2+ 
    #f10 = FFT_ionS(filename=f10, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)#Mass2 and Mass 16 exactly anti phase
    f11 = FFT_ionS(filename=f11, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11)
    #f12 = FFT_ionS(filename=f12, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)
    
    #f29 = FFT_ionS(filename=f29, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11) #MCP 3600
    #f30 = FFT_ionS(filename=f30, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.2,5),ChN=11) #MCP 3800
    #f31 = FFT_ionS(filename=f31, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.2,5),ChN=11)
    #f32 = FFT_ionS(filename=f32, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11) #3-hour measurement
    f33 = FFT_ionS(filename=f33, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.2,5),ChN=11) #repeat 20220601
#
    ##probe 200mW
    #f23 = FFT_ionS(filename=f23, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)
    #f24 = FFT_ionS(filename=f24, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.5,5),ChN=11)
    #f25 = FFT_ionS(filename=f25, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.5,5),ChN=11) #MCP 3200V
    #f26 = FFT_ionS(filename=f26, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,8.5,5),ChN=11) #MCP 3500V 1.3e-7
    #f27 = FFT_ionS(filename=f27, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.5,5),ChN=11)
    #f28 = FFT_ionS(filename=f28, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(100,8.5,5),ChN=11)
    #fall = FFT_ionS(filename=fall, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,8.2,5),ChN=11)

    #same pump probe power
    f38 = FFT_ionS(filename=f38+f38r, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    #f38r = FFT_ionS(filename=f38r, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    #f38rr = FFT_ionS(filename=f38rr, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(150,11.6,5),ChN=11)
    #f39 = FFT_ionS(filename=f39, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(170,10.7,5),ChN=11)#laser drift
    #ftest = FFT_ionS(filename=ftest, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)

    #try
    #f43 = FFT_ionS(filename=f43, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(170,10.7,5),ChN=11)#150-140
    #f44 = FFT_ionS(filename=f44, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(170,10.7,5),ChN=11)#150-160
    #f46 = FFT_ionS(filename=f46, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(170,10.7,5),ChN=11)#150-160
    #f47 = FFT_ionS(filename=f47, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(180,10.7,5),ChN=11)#180-180
    f48 = FFT_ionS(filename=f48, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9.5,5),ChN=11)#200-200
    f49 = FFT_ionS(filename=f49, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,8.5,5),ChN=11)#200-200

    #f50 = FFT_ionS(filename=f50, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9,5),ChN=11)#200-200
    #f51 = FFT_ionS(filename=f51, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(190,9,5),ChN=11)#190-200
    #f52 = FFT_ionS(filename=f52, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(180,9,5),ChN=11)#180-200
    f53 = FFT_ionS(filename=f53, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(180,10,5),ChN=11)#180-180
    #f54 = FFT_ionS(filename=f54, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9,5),ChN=11)#200-200
    #f55 = FFT_ionS(filename=f55, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9,5),ChN=11)#150-150
    #f56 = FFT_ionS(filename=f56, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,9,5),ChN=11)#200-150
####################################################################################
    f0 = FFT_ionS(filename=f0, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(250,9.2,5),ChN=11)#250-250
    f1 = FFT_ionS(filename=f1, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(230,9.7,5),ChN=11)#230-230
    f2 = FFT_ionS(filename=f2, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(200,10.5,5),ChN=11)#200-200
    #ff=[fp1,fp2,fp3,fp4]
    #ff=[f53,f48,f49]#
    ####################################################################################################################
    f00 = FFT_ionS(filename=f00, scanTime=10, sampRate=300, molecule='H2O', intensity=cal_intensity(260,8,5),ChN=11)#260-250

    #ff=[f18+f15+f11+f33]#
    ff=[f38,f53,f48]#measured in June 2022, best data f38,f53,
    for x in ff:
        x.pathFinder()
        print(x.filepath)
        x.read_splitB(isCorrection=True,sumN=5)
        #x.findZeroDelay()
        #x.window(windowSize=150, direction='left')
        #x.rmvExp()
        #x.smooth(windowSize=9)
        #x.show_Spectra(ifsaveT = False)#ifsaveT=True will open a Origin window and export the data. Only used
        #plt.show()
        #x.padding(paddingSize=100000)
        #x.FFTS()
        #x.show_FFT()
        #x.phaseRetrive2()
        
    