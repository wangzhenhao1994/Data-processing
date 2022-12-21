from operator import xor

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
from math import ceil
#from pynufft import NUFFT

my_path = os.path.abspath(__file__)


# font = {
#        'size'   : 35}
#mpl.rc('font', **font)
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Helvetica"],
#     'font.size': 22})
# # for Palatino and other serif fonts use:
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "serif",
#     "font.serif": ["Palatino"],
#     'font.size': 22
# })

SMALL_SIZE = 18
MEDIUM_SIZE = 22
BIGGER_SIZE = 26

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


class FFT_ionS():
    #def __init__(self, filename, massRange=[5, 90], dcRange=2, cal_mass=[45,61], cal_pixel=[415,1141]):
    def __init__(self, filename, scanTime, sampRate, label=None, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[275,469]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[494,583]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[197,279]):
        '''

        '''
        
        self.filename = filename
        self.channelSize = 1024
        self.scanTime = scanTime
        self.label = label
        self.trueScanTime = ceil(self.scanTime/1.6)*1.6
        self.scanLength = self.scanTime*sampRate
        self.peakRange = [-15, 15]  # range of the peak
        self.delay = np.arange(self.scanLength)/self.scanLength*(self.scanTime/self.trueScanTime)*99*2*3.33564*10**-15+40*2*3.33564*10**-15
        self.rebin_delay = None
        self.ifrebin = False
        self.longstage = True
        self.data = 0
        self.massRange = massRange
        self.dcRange = dcRange
        self.windowSize = 100
        self.calculator = Calibration_mass(cal_mass, cal_pixel)
        self.gasBottle = {
            #"44": self.calculator.cal_pixel(44) + self.peakRange,
            #"45": self.calculator.cal_pixel(45) + self.peakRange,
            #"60": self.calculator.cal_pixel(60) + self.peakRange,
            #"CH3OH+": self.calculator.cal_pixel(32) + self.peakRange,
            #"CH2OH+": self.calculator.cal_pixel(31) + self.peakRange,
            #"CHOH+": self.calculator.cal_pixel(30) + self.peakRange,
            #"CHO+": self.calculator.cal_pixel(29) + self.peakRange,
            #"CO+":self.calculator.cal_pixel(28) +self.peakRange,
            #"CH2+":self.calculator.cal_pixel(14) +self.peakRange,
            #"CH3+":self.calculator.cal_pixel(15) +self.peakRange,
            # "N": self.calculator.cal_pixel(14) +self.peakRange,
            # "O": self.calculator.cal_pixel(16)+self.peakRange,
            "H2O+": self.calculator.cal_pixel(18)+self.peakRange,
            # "N2": self.calculator.cal_pixel(28)+self.peakRange,
            # "$O_2$": self.calculator.cal_pixel(32)+self.peakRange,
            # "N": self.calculator.cal_pixel(14) +self .peakRange,
            # "O": self.calculator.cal_pixel(16)+self.peakRange,
            #"H2O": self.calculator.cal_pixel(18)+self.peakRange,
            #"OH-": self.calculator.cal_pixel(17)+self.peakRange,
            # "N2": self.calculator.cal_pixel(28)+self.peakRange,
            # "O2": self.calculator.cal_pixel(32)+self.peakRange,
            # "Ar": self.calculator.cal_pixel(32)+self.peakRange,
            #"\'CH3OH+\'-\'CH2OH+\'": #[0,self.channelSize]
            #np.append(self.calculator.cal_pixel(17) + self.peakRange, self.calculator.cal_pixel(18) + self.peakRange),
            #"\'CH3OH+\'+\'CH2OH+\'": #[0,self.channelSize]
            #np.append(np.append(self.calculator.cal_pixel(17) + self.peakRange, self.calculator.cal_pixel(18) + self.peakRange),0),
        }
        self.spectraBottle = {}
        self.fftS = {}
        self.stftS = {}
        self.spec = None
        self.scanNo = None

    def checkData(self):
        with h5py.File(self.filename, 'r+') as f:
            print(np.array(f['parameters']))
            bigger = True
            i = 0
            while bigger:
                try:
                    a = f['dataD'][i*self.channelSize,0]
                    i+=1
                except ValueError:
                    bigger = False
                    self.scanNo = i
                    f.create_dataset('scanNo', data=i)
                    print('The number of scan is '+str(i)+'!')

    def read_split(self, overNight = False, firstTry = False, sumNo = 20, usefulRange = [0, 5], cu = False):
        with h5py.File(self.filename, 'r+') as f:
            print(f.keys())
            #self.spec=np.array(f['spectrumLog'])
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
                dataSet_created = {'sum0':[0,self.scanNo]}
            for s in dataSet_created.keys():
                e = s in f
                if e and (not firstTry):
                    print('dataset \''+ s +'\' exists!!!')
                    pass
                else:
                    print('Begin to create dataset \''+ s +'\'!')
                    data = np.zeros((self.channelSize, self.scanLength))
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
                            interData = f[s]  # load and overwrite the data
                            interData[...] = data
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
                        print('Begin to overwrite dataset \''+ 'useful' +'\''+'!')
                        interData = f['useful']  # load and overwrite the data
                        interData[...] = self.data

    def FFT(self, t, y):
        n = len(t)
        delta = (max(t) - min(t)) / (n-1)
        k = int(n/2)
        f = np.arange(k) / (n*delta) / 10**12  # frequency unit THz
        Y = abs(sp.fft.fft(y))[:k]
        return np.array([f, Y])

    def mass_spectra(self):

        for gas in self.gasBottle.keys():
            if len(self.gasBottle[gas])==2:
                [pixelMin, pixelMax] = list(map(int, self.gasBottle[gas]))
                self.spectraBottle[gas] = np.sum(
                    self.data[pixelMin:pixelMax, :], 0)
            elif len(self.gasBottle[gas])==4:
                [pixel0, pixel1, pixel2, pixel3] = list(map(int, self.gasBottle[gas]))
                self.spectraBottle[gas] = np.sum(
                    -self.data[pixel0:pixel1, :], 0)+np.sum(
                    self.data[pixel2:pixel3, :], 0)
            else:
                [pixel0, pixel1, pixel2, pixel3, useless] = list(map(int, self.gasBottle[gas]))
                self.spectraBottle[gas] = np.sum(
                    self.data[pixel0:pixel1, :], 0)+np.sum(
                    self.data[pixel2:pixel3, :], 0)

    def window(self, windowSize=0, direction='left'):
        '''
        windowSize is in fs.
        '''
        for gas in self.spectraBottle.keys():
            if 'rebin' not in gas:
                self.spectraBottle[gas], self.delay = self.inter_window(gas, self.delay, windowSize=windowSize, direction=direction)
            else:
                self.spectraBottle[gas], self.rebin_delay = self.inter_window(gas, self.rebin_delay, windowSize=windowSize, direction=direction)
            
    def inter_window(self, gas, delay, windowSize=0, direction='left'):
        '''
        windowSize is in fs.
        '''
        data = self.spectraBottle[gas]
        __len = np.size(data)
        windowSize = int(windowSize/(delay[1]-delay[0])*1e-15)
        print(windowSize)
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
            #self.paddedDelay = np.arange(self.delay[0])
        # if (stepNum_left%2) == 0:
        #    paddingWindow = np.concatenate(( np.zeros(int(stepNum_left/2))+1,
        #    np.zeros(windowSize*2),np.zeros(int(stepNum_left/2))+1,np.zeros(paddingSize)))
        #    #self.paddedDelay = np.arange(self.delay[0])
        # else:
        #    paddingWindow = np.concatenate((np.zeros(paddingSize), np.zeros(int((stepNum_left-1)/2))+1,
        #    np.zeros(windowSize*2),np.zeros(int((stepNum_left+1)/2))+1, np.zeros(paddingSize)))
            # plt.plot(dcShift)
        return window, delay

    def padding(self, paddingSize=0):
        for gas in self.spectraBottle.keys():
            if "rebin" not in gas:
                self.spectraBottle[gas], self.delay = self.inter_padding(gas, self.delay, paddingSize = paddingSize)
            else:
                self.spectraBottle[gas], self.rebin_delay = self.inter_padding(gas, self.rebin_delay, paddingSize = paddingSize)

    def inter_padding(self, gas, delay, paddingSize=0):
        data = None
        delayStep = delay[90]-delay[89]
        delay = np.concatenate((
            np.arange(delay[0]-(paddingSize+1) *
                      delayStep, delay[0], delayStep),
            delay,
            np.arange(delay[-1]+delayStep, delay[-1] +
                      (paddingSize)*delayStep, delayStep)
        ))
        inter_data = self.spectraBottle[gas]
        data = np.concatenate(
            ((np.zeros(paddingSize)), inter_data, (np.zeros(paddingSize))), axis=0)
        delay = delay[:len(data)]
        return data, delay
        

    def interpS(self, interNum=1):
        interpDelay = np.arange(
            self.delay[0], self.delay[-1], (self.delay[1]-self.delay[0])/(interNum+1))
        for gas in self.spectraBottle.keys():
            iS = np.interp(interpDelay, self.delay, self.spectraBottle[gas])
            self.spectraBottle[gas] = iS
            # print(np.shape(iS))
        self.delay = interpDelay
        # print(np.sie(self.delay))

    def FFTS(self):
        for gas in self.spectraBottle.keys():
            if self.ifrebin:
                if "rebin" not in gas:
                    continue
                else:
                    delay = self.rebin_delay
                    self.fftS[gas] = self.FFT(
                        delay, self.spectraBottle[gas][self.windowSize::])
            else:
                if "rebin" not in gas:
                    delay = self.delay
                    self.fftS[gas] = self.FFT(
                        delay, self.spectraBottle[gas][self.windowSize::])
                else:
                    continue

    # def filterS(self, band):
     #   for gas in self.spectraBottle.keys():
      #      N, wn = sps.buttord(1e11, )
       #     fil = sps.butter(10, band, btype='band', fs=1000, output='sos')
        #    self.spectraBottle[gas] = sps.sosfilt(fil, self.spectraBottle[gas])

    def STFTS(self, windowsize=0, ratio=10):
        '''
        windowsize is window size in fs
        '''
        for gas in self.spectraBottle.keys():
            #if 'rebin' not in gas or 'filter' not in gas:
            if 'rebin' not in gas:
                continue
            f, t, self.stftS[gas] = sps.stft(self.spectraBottle[gas], fs=1/(
                self.rebin_delay[2]-self.rebin_delay[1]), noverlap=windowsize-2, nperseg=windowsize, nfft=windowsize*3)
            #self.stftS[gas] = self.stftS[gas][:,::-1]*10
            self.stftS[gas] = self.stftS[gas]*10
            vmax=abs(self.stftS[gas][0:int(len(f)/10)]).max()
            vmin=abs(self.stftS[gas]).min()
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            levels = np.arange(vmin,vmax,(vmax-vmin)/500)
            plt.figure(figsize=(8.0, 5.0))
            plt.contourf(t*10**15, f/10**12*33.35641,
                           np.abs(self.stftS[gas]), levels=levels, cmap='jet', norm=norm)
            #plt.xlim([100, self.delay[-1]])
            plt.ylim([000, 4500])
            plt.ylabel('Frequency [cm-1]')
            plt.xlabel('Time [fs]')
            #plt.clim(0,abs(self.stftS[gas][:int(len(f)/15)]).max())
            plt.colorbar(ticks=levels[::100])
            d = '_55_90'
            #plt.title(gas+d)
            plt.tight_layout()
            plt.savefig("20211103_sft_18"+d,dpi=720, bbox_inches='tight', pad_inches=0.2)
            plt.show()

    def rmvExp(self):
        for gas in self.spectraBottle.keys():
            #y = self.spectraBottle[gas]
            #k, b = np.polyfit(f.delay, np.log(y), 1)#, w=np.sqrt(y))
            #self.spectraBottle[gas] = y- np.exp(b)*np.exp(f.delay*k)
            self.spectraBottle[gas]=polynomial(self.spectraBottle[gas], order=2, plot=False)

    def runingAverage(self, n=5):
        def runAve(x, N): return np.convolve(x, np.ones(N)/N, mode='valid')
        for gas in self.spectraBottle.keys():
            self.spectraBottle[gas] = runAve(self.spectraBottle[gas], n)
            new_size = len(self.spectraBottle[gas])
        self.delay = self.delay[:new_size]

    def smooth(self, windowSize=100, order=3):
        for gas in self.spectraBottle.keys():    
            self.spectraBottle[gas]=sps.savgol_filter(self.spectraBottle[gas], windowSize, order) # window size 51, polynomial order 3
    
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
        for gas in list(self.spectraBottle.keys()):
            self.spectraBottle['rebin_'+gas] = self.rebin_factor(self.spectraBottle[gas], factor)
        self.rebin_delay = self.rebin_factor(self.delay, factor)

    def show_FFT(self):
        #plt.figure()
        #unused_f, ref = self.fftS['CH3OH+']
        for gas in self.fftS.keys():
            [f, Y] = self.fftS[gas]
            plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/np.sum(Y[self.dcRange:]), label=self.label)#float("{:.2f}".format(self.trueScanTime)))
            plt.xlabel('Frequency/cm-1')
            plt.ylabel('a.u.')
            plt.legend(gas)
            #plt.ylim([0,np.max(Y[self.dcRange:1000])*3/2])
            plt.xlim([0,20000])
        plt.legend()
        plt.tight_layout()
        plt.title('Pump intensity scan_Area normalization')
        #plt.savefig("fft_31"+d,dpi=720, bbox_inches='tight', pad_inches=0.2)

        # plt.figure()
        # for gas in ['N', 'O']:
        #    [f, Y] = self.fftS[gas]
        #    plt.plot(f[self.dcRange:], Y[self.dcRange:], label=gas)
        #    plt.xlabel('THz')
        #    plt.ylabel('a.b.')
        #    plt.legend(gas)
        # plt.legend()
        #plt.show()

    def show_Spectra(self, shift=0):
        #plt.figure()
        for gas in self.spectraBottle.keys():
            #if 'rebin' not in gas or 'filter' not in gas:
            #if 'filter' not in gas:
            #    continue
            #elif 'rebin' not in gas:
            delay = self.delay
            #else:
            #    delay = self.rebin_delay
            plt.plot(delay*10**15+shift,
                     self.spectraBottle[gas])
            #plt.plot(self.spectraBottle[gas]/np.sum(self.spectraBottle[gas]), label=gas)
            plt.xlabel("Delay/fs")
            plt.ylabel('a.u.')
            #plt.xlim([200,600])
            #plt.legend(gas[:])
        #plt.savefig("spectra_31"+d,dpi=720, bbox_inches='tight', pad_inches=0.2)

    def dataProcess(self):
        # self.mass_spectra()
        self.FFTS()

    # def filter(self, lowcut, highcut, fs, order=5):
    #     #b, a = self.butter_bandpass(lowcut, highcut, fs, order=order)
    #     #w, h = sps.freqz(b, a, worN=2000)

    #     sos = self.butter_bandpass(lowcut, highcut, fs, order=order)
    #     w, h = sps.sosfreqz(sos, worN=2000)

    def useFilter(self, lowcut, highcut):
        for gas in list(self.spectraBottle.keys()):
            # if 'rebin' not in gas:
            fs = 1/(self.delay[90]-self.delay[89])
            # else:
            #     fs = 1/(self.rebin_delay[90]-self.rebin_delay[89])
            self.spectraBottle['filter_'+gas] = self.butter_bandpass_filter(self.spectraBottle[gas], lowcut, highcut, fs)

    def butter_bandpass_filter(self, data, lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = sps.butter(order, [low, high], analog=False, btype='band', output='sos')
        y = sps.sosfiltfilt(sos, data)
        #w, h = sps.sosfreqz(sos, worN=2000)
        #plt.plot(w, h)
        #plt.show()
        return y

class FFTOSC(FFT_ionS):
    def __init__(self, delay, data, massRange=[31, 32], dcRange=2, label=None):
        super().__init__(delay, data, massRange=massRange, dcRange=dcRange, label=label)

    def mass_spectra(self):
        self.spectraBottle['CH3OH+'] = self.data


if __name__ == '__main__':
    
    dataPath = pl.PureWindowsPath(
        r'D:\DataProcessing\2022-07\2022-07-27')
    
    #datafile =  os.path.join(dataPath, r'scan_tof_2021-12-19-18-14-47.hdf5') #70-90 time:10
    ##datafile2 =  os.path.join(dataPath, r'scan_tof_2021-12-19-19-38-40.hdf5') #75-90 time:10
    #datafile3 =  os.path.join(dataPath, r'scan_tof_2021-12-19-22-50-06.hdf5') #80-90 time:10
    #datafile4 =  os.path.join(dataPath, r'scan_tof_2021-12-20-09-44-01.hdf5') #90-90 time:10
    #datafile5 =  os.path.join(dataPath, r'scan_tof_2021-12-20-11-17-06.hdf5') #100-90 time:10
    #datafile5_2 =  os.path.join(dataPath, r'scan_tof_2021-12-20-12-56-02.hdf5') #100-90 time:16
    
    datafile6 =  os.path.join(dataPath, r'scan_tof_2022-07-27-12-17-05.hdf5') #110-90 time:10
    

    
    
    ##f = FFT_ionS(filename=datafile, scanTime=200, sampRate=1000, dcRange=2, massRange=[5,50])#longstage
    #f = FFT_ionS(filename=[datafile], scanTime=10, label="70pp90pr", sampRate=1000, dcRange=2, massRange=[5,50])
    ##f2 = FFT_ionS(filename=[datafile2], scanTime=10, sampRate=1000, dcRange=2, massRange=[5,50])
    #f3 = FFT_ionS(filename=[datafile3], scanTime=10, label="80pp90pr", sampRate=1000, dcRange=2, massRange=[5,50])
    #f4 = FFT_ionS(filename=[datafile4], scanTime=10, label="90pp90pr", sampRate=1000, dcRange=2, massRange=[5,50])
    #f5 = FFT_ionS(filename=[datafile5], scanTime=10, label="100pp90pr", sampRate=1000, dcRange=2, massRange=[5,50])
    #f5_2 = FFT_ionS(filename=[datafile5_2], scanTime=16, sampRate=1000, dcRange=2, massRange=[5,50])
    f6 = FFT_ionS(filename=datafile6, scanTime=10, label="110pp90pr", sampRate=1000, dcRange=2, massRange=[5,50])
    

    ff=[f6]
    #ff=[f5_2]

    for x in ff:
        x.spectraBottle['H2O+'] = 0
        x.read_split()
        for file in x.filename:
            with h5py.File(file, 'r+') as f:
                print(f.keys())
                #print(f['parameters'].value)
                print(f['dataD'].shape)
                x.spectraBottle['H2O+'] = x.spectraBottle['H2O+']+np.sum(np.array(f['dataD'], dtype=float),axis=0)
                x.spectraBottle['H2O+'] = x.spectraBottle['H2O+']
                print(x.spectraBottle['H2O+'].shape)
        
        x.window(windowSize=50, direction='left')
        #pnf.naff(signal, 500, 1, 0 , False)
        #x.smooth(windowSize=21)
        #x.rmvExp()
        #x.rebinS(factor=10)
        #x.padding(paddingSize=150000)
        #x.STFTS(windowsize=int(150/600*2000))
        #x.dataProcess()
        #x.show_Spectra()
        #plt.show()
        #
        #x.show_FFT()
    #plt.show()


        
        
    