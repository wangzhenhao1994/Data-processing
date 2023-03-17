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
import re
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
    def __init__(self, filename, scanTime, sampRate, massRange=[5, 50], molecule='CO2', intensity = 0, dcRange=5, ChN=11, cal_mass=[2,16], cal_pixel=[3082,10216], direction = 'left'):
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
        self.direction = direction
        try:
            self.intensity = 'pu'+str('%.1E' % Decimal(intensity[0]))+'pr'+str('%.1E' % Decimal(intensity[1]))
        except TypeError:
            self.intensity = str('%.1E' % Decimal(intensity))
        self.stage = 'piezo'

        self.rootPath= pl.PureWindowsPath(r'D:\SF_FFT')
        self.savePath = os.path.join(self.rootPath, r'processedData', r'inter', self.molecule, str(self.intensity)+'_'+self.molecule)


    def pathFinder(self):
        for fn in self.filename:
            interpath = pl.PureWindowsPath(fn[9:16], fn[9:19])
            self.filepath = self.filepath + [pl.PureWindowsPath(self.rootPath, r'OriginalData', interpath, fn + r'.hdf5')]
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
            self.delayB = np.arange(self.scanLengthB)/self.scanLengthB*100*2*2*3.33564*10**-15*3642/3647
            if isCorrection:
                with h5py.File(fp, 'r+') as f:
                    print(f.keys())
                    try:
                        data = np.array(f['dataB'], dtype=float)
                    except OSError:
                        self.num=0
                        continue
                    print(data.shape)
                    m = int(data.shape[0]/self.ChN/sumN)
                    print('m is '+ str(m)+' !')
                    for i in [0,2,4,6,8,10]:
                        self.interSpectraBottleB['Ch'+str(i)] = data[i::ChN][-int(m*sumN):].reshape(m, sumN, data.shape[1]).sum(axis=1)
                    self.interinterSpectraBottleB = {}
                    num = int(m/10)
                    self.num = num
                    if num ==0:
                        continue
                    for gas in self.interSpectraBottleB.keys():
                        self.interinterSpectraBottleB[gas] = {}
                        self.spectraBottleB[gas] = np.zeros((num, 13000))
                    for i in range(num):
                        for gas in self.interSpectraBottleB.keys():
                            _inter = self.interSpectraBottleB[gas][i*10:(i+1)*10]
                            if self.direction == 'right':
                                self.interinterSpectraBottleB[gas][str(i)] = np.flip(_inter,axis=1)
                            else:
                                self.interinterSpectraBottleB[gas][str(i)] = _inter
            
            self.delayCorrection()
            if not os.path.exists(self.savePath):
                os.mkdir(self.savePath)
            print(self.interSpectraBottleB['Ch0'].shape)
            save_obj(self.interSpectraBottleB, pl.PureWindowsPath(self.savePath, os.path.split(fp)[1].replace(r'.hdf5','')+r'.pkl'))

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

    def calDrift(self, _cRange):
        plt.clf()
        _iS = np.array(np.zeros(self.interSpectraBottleB['Ch0'].shape))
        sumSpec=_iS
        for gas in self.interSpectraBottleB.keys():
            sumSpec=sumSpec+self.interSpectraBottleB[gas]
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            print(i)
            #_iS[i]=self.butter_bandpass_filter(self.interSpectraBottleB[gas][i], (1389.15-10)/33.35641*1e12, (1389.15+10)/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
            _iS[i]=self.butter_bandpass_filter(sumSpec[i], 500/33.35641*1e12, 2000/33.35641*1e12, 1/(self.delayB[1]-self.delayB[0]))
        iS = sp.interpolate.RectBivariateSpline(range(_iS.shape[0]), self.delayB*1e15, _iS)
        _delayRange = np.linspace(start=_cRange[0],stop=_cRange[1], num=20000)
        indexMax = []
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange)
            indexMax = indexMax + [np.argmax(np.abs(_inter))]
            #plt.plot(_delayRange,_inter)
        
        print(indexMax)
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        #plt.show()
        _ref = sum(indexMax[int(self.interSpectraBottleB[gas].shape[0]/2)-5:int(self.interSpectraBottleB[gas].shape[0]/2)+5])/10
        _shift = (np.array(indexMax)-_ref)*(_cRange[1]-_cRange[0])/20000
        for i in range(self.interSpectraBottleB[gas].shape[0]):
            _inter=iS.ev(i,_delayRange+_shift[i])
            #plt.plot(_delayRange,_inter)
        #plt.xlabel('Delay (fs)')
        #plt.ylabel('a.u.')
        #plt.show()
        return _shift

    def delayCorrection(self, _cRange = [0,150]):
        for k in range(self.num):
            for gas in self.interSpectraBottleB.keys():
                self.interSpectraBottleB[gas] = self.interinterSpectraBottleB[gas][str(k)]
            xxx = 0
            while True:
                try:
                    print(_cRange)
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
                if gas == 'Ch8':
                    print(sum(inter))
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
                     self.spectraBottleB[gas][0], label=gas)#/np.amax(self.spectraBottleB[gas])
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
        plt.savefig(pl.PureWindowsPath(self.savePath, self.filename[0]+r'.png'),dpi=720, bbox_inches='tight', pad_inches=0.2)
        #plt.show()

if __name__ == '__main__':

    directory = os.path.abspath(r'D:\SF_FFT\OriginalData')
    ff=[]
    for root, dirs, files in os.walk(directory):
        #print(root)
        path = root.split(os.sep)
        for file in files:
            if '.txt' in file:
                with open(os.path.join(root, file), 'r') as f:
                    a = f.read()
                    if r'Ch8 -> Mass 44 CO2+' in a:
                        fname = file.replace('.txt','')
                        interIn=list(map(int, re.findall(r'\d+', a[0:45])))
                        inter = list(map(float, re.findall(r"[-+]?(?:\d*\.\d+|\d+)", a[45:70])))
                        if interIn.__len__()==2 and inter.__len__()==2:
                            intensity=[cal_intensity(int(interIn[0]),float(inter[0]),5),cal_intensity(int(interIn[1]),float(inter[1]),5)]
                        elif interIn.__len__()==1 and inter.__len__()==2:
                            intensity=[cal_intensity(int(interIn[0]),float(inter[0]),5),cal_intensity(int(interIn[0]),float(inter[1]),5)]
                        elif interIn.__len__()==2 and inter.__len__()==1:
                            intensity=[cal_intensity(int(interIn[0]),float(inter[0]),5),cal_intensity(int(interIn[1]),float(inter[0]),5)]
                        elif interIn.__len__()==1 and inter.__len__()==1:
                            intensity=[cal_intensity(int(interIn[0]),float(inter[0]),5),cal_intensity(int(interIn[0]),float(inter[0]),5)]
                        else:
                            print('Somethig is wrong!\n')

                        if r'overlap at' in a:
                            direction='right'
                            intensity.reverse()
                        else:
                            direction='left'
                        ff=ff+[FFT_ionS(filename=[fname], scanTime=10, sampRate=300, molecule='CO2', intensity=intensity,ChN=11,direction=direction)]
    
        #CO2
    #f3 = [r'scan_tof_2022-07-27-12-17-05'] #230-230
    #f4 = [r'scan_tof_2022-07-27-13-30-56'] #150-150
    #f5 = [r'scan_tof_2022-07-27-16-38-24'] #120-180
    #f6 = [r'scan_tof_2022-07-27-16-52-35'] #260-120
    #f7 = [r'scan_tof_2022-07-28-15-44-21'] #250-70 r'scan_tof_2022-07-27-18-33-55',
    #f8 = [r'scan_tof_2022-07-28-11-10-28'] #210-70
    #f9 = [r'scan_tof_2022-07-28-12-33-51'] #170-70
    #f10 = [r'scan_tof_2022-07-28-14-07-58'] #130-70
#
    #f11 = [r'scan_tof_2022-07-29-13-10-14'] #70-250
    #f12 = [r'scan_tof_2022-07-29-14-48-30'] #70-210
    #f13 = [r'scan_tof_2022-07-29-16-22-58'] #70-170
#
#
    #f14 = [r'scan_tof_2022-07-28-17-18-03'] #170-120
    #f15 = [r'scan_tof_2022-07-28-18-39-13'] #170-150
    #f16 = [r'scan_tof_2022-07-28-20-04-15'] #170-180
    #f17 = [r'scan_tof_2022-07-29-09-52-06'] #170-210
    #f18 = [r'scan_tof_2022-07-29-11-23-24'] #170-250
#
    #f19 = [r'scan_tof_2022-07-29-17-42-30'] #250-270
#
    #f7 = FFT_ionS(filename=f7, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(250,9,5),cal_intensity(70,19.5,5)],ChN=11)
    #f8 = FFT_ionS(filename=f8, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(210,10,5),cal_intensity(70,19.5,5)],ChN=11)
    #f9 = FFT_ionS(filename=f9, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(170,12,5),cal_intensity(70,19.5,5)],ChN=11)
    #f10= FFT_ionS(filename=f10, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(130,14.5,5),cal_intensity(70,19.5,5)],ChN=11)
    #f11= FFT_ionS(filename=f11, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(70,19.5,5),cal_intensity(250,9,5)],ChN=11)
    #f12= FFT_ionS(filename=f12, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(70,19.5,5),cal_intensity(210,10,5)],ChN=11)
    #f13= FFT_ionS(filename=f13, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(70,19.5,5),cal_intensity(170,12,5)],ChN=11)
#
    #f14= FFT_ionS(filename=f14, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(170,12,5),cal_intensity(120,13.5,5)],ChN=11)
    #f15= FFT_ionS(filename=f15, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(170,12,5),cal_intensity(150,12,5)],ChN=11)
    #f16= FFT_ionS(filename=f16, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(170,12,5),cal_intensity(180,10.5,5)],ChN=11)
    #f17= FFT_ionS(filename=f17, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(170,12,5),cal_intensity(210,10,5)],ChN=11)
    #f18= FFT_ionS(filename=f18, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(170,12,5),cal_intensity(250,9,5)],ChN=11)
    #f19= FFT_ionS(filename=f19, scanTime=10, sampRate=300, molecule='CO2', intensity=[cal_intensity(250,8.4,5),cal_intensity(270,8.4,5)],ChN=11)

    #ff=[fp1,fp2,fp3,fp4]
    #ff=[f19]#,f11,f12,f13,f14,f15,f16,f17,f18,f19
    for x in ff:
        x.pathFinder()
        print(x.filepath)
        
        x.read_splitB(isCorrection=True,sumN=5)
        if x.num==0:
            continue
        #x.findZeroDelay()
        #x.window(windowSize=150, direction='left')
        #x.rmvExp()
        #x.smooth(windowSize=9)
        x.show_Spectra()
        #plt.show()
        #x.padding(paddingSize=100000)
        #x.FFTS()
        #x.show_FFT()
        #x.phaseRetrive2()
        
    