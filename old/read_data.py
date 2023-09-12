import sys

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

    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[50,62]):
    def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[1,18], cal_pixel=[730,9399]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[197,279]):
        '''

        '''
        
        self.filename = filename
        self.channelSize = 12032
        self.scanLength = 3320#72000
        self.peakRange = [-15, 15]  # range of the peak
        self.delay = np.arange(self.scanLength)/self.scanLength*(116.7/119)*(119.5/120*99)*2*3.33564*10**-15
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
            "H2O++": self.calculator.cal_pixel(9)+self.peakRange,
            # "N2": self.calculator.cal_pixel(28)+self.peakRange,
            # "$O_2$": self.calculator.cal_pixel(32)+self.peakRange,
            # "N": self.calculator.cal_pixel(14) +self .peakRange,
            # "O": self.calculator.cal_pixel(16)+self.peakRange,
            #"H2O": self.calculator.cal_pixel(18)+self.peakRange,
            #"OH+": self.calculator.cal_pixel(17)+self.peakRange,
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
            #print(np.array(f['ppNum']))
            print(f['dataD'].shape)
            bigger = True
            i = 0
            while bigger:
                try:
                    a = f['dataD'][i*self.channelSize,0]
                    i+=1
                except IndexError:
                    bigger = False
                    self.scanNo = i
                    f.create_dataset('scanNo', data=i)
                    print('The number of scan is '+str(i)+'!')

    def read_split(self, overNight = False, firstTry = False, sumNo = 20, usefulRange = [0, 5], cu = False):
        with h5py.File(self.filename, 'r+') as f:
            #plt.plot(np.sum(np.array(f['dataD'][164*self.channelSize:(164+1)*self.channelSize]),0))
            #print('the pump-only count is '+str(np.sum(np.array(f['dataD'][164*self.channelSize:(164+1)*self.channelSize]))/116.7),'\nthe probe-only count is '+str(np.sum(np.array(f['dataD'][165*self.channelSize:(165+1)*self.channelSize]))/116.7) )
            #plt.show()
            print(f.keys())
            #print(np.array(f['ppNum']))
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
                    except ValueError:
                        print('Dataset \''+ s +'\''+' already exists!')
                        try:
                            interData = f[s]  # load and overwrite the data
                            interData[...] = data
                            print('Dataset \'' + s + '\''+' is overwrited!')
                        except ValueError:
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
                except ValueError:
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
            vmax=abs(self.stftS[gas][int(len(f)/40):int(len(f)/15)]).max()
            vmin=abs(self.stftS[gas]).min()
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            levels = np.arange(vmin,vmax,(vmax-vmin)/500)
            plt.figure(figsize=(8.0, 5.0))
            plt.contourf(t*10**15+633, f/10**12*33.35641,
                           np.abs(self.stftS[gas]), levels=levels, cmap='jet', norm=norm)
            #plt.xlim([100, self.delay[-1]])
            plt.ylim([1000, 4500])
            plt.ylabel('Frequency [cm-1]')
            plt.xlabel('Time [fs]')
            plt.clim(0,abs(self.stftS[gas][:int(len(f)/15)]).max())
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
            self.spectraBottle[gas]=polynomial(self.spectraBottle[gas], order=6, plot=False)

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
            plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:], label=gas)
            plt.xlabel('Frequency/cm-1')
            plt.ylabel('a.u.')
            plt.legend(gas)
            #plt.ylim([0,np.max(Y[self.dcRange:1000])*3/2])
            plt.xlim([0,22000])
        plt.legend()
        plt.tight_layout()
        d='_100_120'
        #plt.title('FFT[100:600]'+d)
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
                     self.spectraBottle[gas], label=gas)
            #plt.plot(self.spectraBottle[gas]/np.sum(self.spectraBottle[gas]), label=gas)
            plt.xlabel("Delay/fs")
            plt.ylabel('a.u.')
            #plt.xlim([200,600])
            plt.legend(gas)
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
    
    if 'linux' in sys.platform:
        dataPath = r'/mnt/c/Users/user/Desktop/Data_newTOF'
        dataPath = pl.PurePath(dataPath)
    else:
        dataPath = r'D:\SF_FFT\OriginalData'
        dataPath = pl.PureWindowsPath(dataPath)
    yearAndMonth = r'2022-07'
    
    
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-11', r'scan_tof_2021-10-11-23-23-40', False, 20, [5,11], True]#120-120 80 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-16', r'scan_tof_2021-10-16-19-47-59', False, 85, [0,1], False]#120-120 80 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-14', r'scan_tof_2021-10-14-21-00-11', False, 36, [0,1], False]#100-110 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-14', r'scan_tof_2021-10-14-19-17-24', False, 36, [0,1], False]#90-110 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-14', r'scan_tof_2021-10-14-17-50-39', False, 36, [0,1], False]#80-110 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-14', r'scan_tof_2021-10-14-16-08-58', False, 36, [0,1], False]#70-110 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-14', r'scan_tof_2021-10-14-14-40-33', False, 36, [0,1], False]#60-110 35 Ar
    
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-15', r'scan_tof_2021-10-15-11-22-53', False, 36, [0,1], False]#120-90 35 Ar # two CEP, try again #no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-15', r'scan_tof_2021-10-15-14-57-49', False, 36, [0,1], False]#110-90 35 Ar# two CEP, try again #no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-15', r'scan_tof_2021-10-15-13-00-34', False, 36, [0,1], False]#110-90 35 Ar# two CEP, try again #no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-14', r'scan_tof_2021-10-14-13-05-18', False, 36, [0,1], False]#100-90 35 Ar #no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-14', r'scan_tof_2021-10-14-11-20-21', False, 36, [0,1], False]#90-90 35 Ar #no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-15', r'scan_tof_2021-10-15-18-41-16', False, 36, [0,1], False]#85-90 35 Ar #no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-14', r'scan_tof_2021-10-14-09-48-02', False, 36, [0,1], False]#80-90 35 Ar #spectra has blue drift,try again later #no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-15', r'scan_tof_2021-10-15-16-56-57', False, 36, [0,1], False]#75-90 35 Ar, two peaks
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-13', r'scan_tof_2021-10-13-19-52-07', False, 36, [0,1], False]#70-90 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-13', r'scan_tof_2021-10-13-18-12-42', False, 36, [0,1], False]#60-90 35 Ar

    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-13', r'scan_tof_2021-10-13-16-36-09', False, 36, [0,1], False]#100-100 35 Ar try again later, saturation.
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-12', r'scan_tof_2021-10-13-12-49-40', False, 36, [0,1], False]#90-100 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-15', r'scan_tof_2021-10-15-20-58-30', False, 36, [0,1], False]#85-100 35 Ar, no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-13', r'scan_tof_2021-10-13-14-53-44', False, 36, [0,1], False]#80-100 35 Ar, strong signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-12', r'scan_tof_2021-10-12-20-17-16', False, 36, [0,1], False]#80-100 35 Ar, strong signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-16', r'scan_tof_2021-10-16-13-15-03', False, 36, [0,1], False]#80-100 35 Ar, no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-16', r'scan_tof_2021-10-16-15-19-13', True, 36, [0,1], True]#80-100 35 Ar, little signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-16', r'scan_tof_2021-10-16-16-57-38', False, 36, [0,1], False]#90-110 35 Ar, no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-16', r'scan_tof_2021-10-16-11-29-18', False, 36, [0,1], False]#75-100 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-13', r'scan_tof_2021-10-13-11-05-44', False, 36, [0,1], False]#70-100 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-13', r'scan_tof_2021-10-13-09-20-24', False, 36, [0,1], False]#60-100 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-12', r'scan_tof_2021-10-12-21-58-45', False, 15, [0,4], False]#60-100 250 Ar

    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-13', r'scan_tof_2021-10-13-21-40-53', False, 15, [0,5], False]#60-80 250 Ar #no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-15', r'scan_tof_2021-10-15-09-39-43', False, 36, [0,1], False]#110-70 35 Ar #no signal # two CEP, try again
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-14', r'scan_tof_2021-10-14-22-57-02', False, 15, [0,4], False]#110-60 150 Ar #no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-12', r'scan_tof_2021-10-12-16-05-08', False, 30, [0,1], False]

    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-03', r'scan_tof_2021-11-03-19-00-50', False, 63, [0,1], False]#55-90 60, strong signal 3600
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-03', r'scan_tof_2021-11-03-22-21-04', False, 10, [0,4], False]#65-90 90, no signal, low count due to pressure decrease
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-04', r'scan_tof_2021-11-04-09-52-26', False, 41, [0,1], False]#75-90 40, 1600cm-1?
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-04', r'scan_tof_2021-11-04-15-10-58', False, 41, [0,1], False]#75-90, noise due to the door
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-04', r'scan_tof_2021-11-04-16-52-55', False, 41, [0,1], False]#75-90
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-04', r'scan_tof_2021-11-04-18-55-52', False, 41, [0,1], False]#75-90
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-04', r'scan_tof_2021-11-04-22-30-00', False, 32, [0,3], False]#85-90 90
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-05', r'scan_tof_2021-11-05-10-16-32', False, 34, [0,1], False]#95-90 35 nothjing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-05', r'scan_tof_2021-11-05-11-52-30', False, 34, [0,1], False]#105-90 35 1600?
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-05', r'scan_tof_2021-11-05-13-49-30', False, 34, [0,1], False]#115-90 35 medium signal 3600
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-05', r'scan_tof_2021-11-05-15-18-14', False, 34, [0,1], False]#125-90 35 medium signal 3550
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-05', r'scan_tof_2021-11-05-16-59-00', False, 34, [0,1], False]#135-90 35 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-05', r'scan_tof_2021-11-05-18-40-42', False, 34, [0,1], False]#145-90 35 1550?
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-05', r'scan_tof_2021-11-05-20-41-03', False, 34, [0,1], False]#145-100 35 max pump count, nothing

    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-05', r'scan_tof_2021-11-05-22-53-43', False, 20, [0,10], False]#110-145 200 max pump count, right
    
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-07', r'scan_tof_2021-11-07-11-12-25', False, 63, [0,1], False]#125-100 60 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-07', r'scan_tof_2021-11-07-14-15-37', False, 34, [0,1], False]#125-90 60 nothing
    
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-13', r'scan_tof_2021-11-13-18-30-49', False, 30, [0,1], False]#110-90 30 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-13', r'scan_tof_2021-11-13-20-27-00', False, 34, [0,1], False]#100-90 35 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-13', r'scan_tof_2021-11-13-22-15-02', False, 20, [0,3], False]#120-90 250 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-14', r'scan_tof_2021-11-14-15-03-49', False, 34, [0,1], False]#70-90 35 strong signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-14', r'scan_tof_2021-11-14-16-18-22', False, 34, [0,1], False]#70-90 35 400 fs further, strong signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-14', r'scan_tof_2021-11-14-17-33-29', False, 23, [0,1], False]#70-90 35 800 fs further, basically no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-14', r'scan_tof_2021-11-14-19-11-19', False, 34, [0,1], False]#130-90 35 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-14', r'scan_tof_2021-11-14-20-50-22', False, 34, [0,1], False]#120-120 35 weak signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-14', r'scan_tof_2021-11-14-22-22-22', False, 34, [0,1], False]#130-120 35 weak signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-14', r'scan_tof_2021-11-14-23-56-59', False, 20, [0,1], True]#140-120 250 strong signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-15', r'scan_tof_2021-11-15-10-37-01', False, 34, [0,1], False]#140-140 35 very weak signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-15', r'scan_tof_2021-11-15-11-59-41', False, 34, [0,1], False]#140-130 35 very weak signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-15', r'scan_tof_2021-11-15-13-41-08', False, 34, [0,1], False]#140-150 35 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-15', r'scan_tof_2021-11-15-15-10-56', False, 34, [0,1], False]#140-110 35 strong signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-15', r'scan_tof_2021-11-15-18-32-35', False, 21, [0,1], False]#140-115 30 test openloop, failure
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-15', r'scan_tof_2021-11-15-19-27-05', False, 34, [0,1], False]#140-115 35 medium signal, overlap measured
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-15', r'scan_tof_2021-11-15-22-00-28', False, 20, [0,3], True]#140-115 250 weak signal, 20um away from overlap
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-16', r'scan_tof_2021-11-16-11-20-03', False, 34, [0,1], False]#140-120 35 medium signal, 1500 cm-1,  20um away from overlap, increase laser intensity, better broadening
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-16', r'scan_tof_2021-11-16-12-45-27', False, 34, [0,1], False]#140-120 35 medium signal, nothing, 20um away from overlap
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-16', r'scan_tof_2021-11-16-14-57-31', False, 34, [0,1], False]#110-120 35  weak signal,max pump count 20um away from overlap
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-16', r'scan_tof_2021-11-16-16-17-03', False, 34, [0,1], False]#110-120 35  no signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-16', r'scan_tof_2021-11-16-18-05-40', False, 34, [0,1], False]#110-120 35  no signal,max pump count 20um away from overlap
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-16', r'scan_tof_2021-11-16-20-37-22', False, 34, [0,1], False]#100-120 35  very strong signal,max pump count 20um away from overlap,adjust the noth filter again
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-16', r'scan_tof_2021-11-16-22-39-01', False, 29, [0,1], False]#105-120 30  nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-17', r'scan_tof_2021-11-17-00-04-56', False, 20, [0,10], True]#100-120 250 strong signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-17', r'scan_tof_2021-11-17-11-15-41', False, 34, [0,1], False]#95-120 35 strong signal, 2 peaks, 40um away from overlap
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-17', r'scan_tof_2021-11-17-13-09-19', False, 34, [0,1], False]#95-120 35 nothing, spectra drift, evacuate the fiber tube afterwards
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-17', r'scan_tof_2021-11-17-15-08-13', False, 34, [0,1], False]#95-120 35 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-17', r'scan_tof_2021-11-17-16-45-12', False, 34, [0,1], False]#102-120 35 weird signal, basically nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-17', r'scan_tof_2021-11-17-18-52-25', False, 34, [0,1], False]#99-120 35 nothing

    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-18', r'scan_tof_2021-11-18-12-23-43', False, 34, [0,1], False]#110-120 35 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-18', r'scan_tof_2021-11-18-14-10-08', False, 34, [0,1], False]#100-120 35 nothing
    
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-19', r'scan_tof_2021-11-19-14-09-04', False, 29, [0,1], False]#100-80 35 nothing 132fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-19', r'scan_tof_2021-11-19-15-44-36', False, 29, [0,1], False]#90-80 35 strong signal 132fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-19', r'scan_tof_2021-11-19-17-06-15', False, 34, [0,1], False]#85-80 35 strong signal 132fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-19', r'scan_tof_2021-11-19-18-39-01', False, 34, [0,1], False]#80-80 35 strong signal 132fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-19', r'scan_tof_2021-11-19-20-12-57', False, 34, [0,1], False]#87.5-80 35 strong signal 132fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-19', r'scan_tof_2021-11-19-22-05-44', False, 20, [0,15], False]#87.5-80 250 strong signal 264fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-20', r'scan_tof_2021-11-20-13-43-01', False, 34, [0,1], False]#90-80 35 weak signal 132fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-20', r'scan_tof_2021-11-20-15-19-22', False, 34, [0,1], False]#90-80 35 nothing 132fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-20', r'scan_tof_2021-11-20-16-51-17', False, 34, [0,1], False]#85-80 35 strong signal 132fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-20', r'scan_tof_2021-11-20-18-12-28', False, 34, [0,1], False]#85-80 35 strong signal 633fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-20', r'scan_tof_2021-11-20-20-15-37', False, 34, [0,1], False]#82.5-80 35 strong signal 132fs
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-20', r'scan_tof_2021-11-20-22-04-00', False, 30, [0,5], True]#85-70 300 strong signal 633fs

    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-21', r'scan_tof_2021-11-21-15-04-43', False, 34, [0,1], False]#60-80 35 WEAK signal 132fs Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-21', r'scan_tof_2021-11-21-16-50-58', False, 34, [0,1], False]#60-80 35 weaki signal 132fs Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-11-21', r'scan_tof_2021-11-21-18-22-14', False, 34, [0,1], False]#50-80 35 weaki signal 132fs Ar

    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-12-08', r'scan_tof_2021-12-08-11-51-27', False, 34, [0,1], False]#85-80 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-12-08', r'scan_tof_2021-12-08-13-39-54', False, 29, [0,1], False]#80-80 nothing
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-12-08', r'scan_tof_2021-12-08-15-05-51', False, 25, [0,1], False]#90-80 strong signal
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-12-08', r'scan_tof_2021-12-08-16-28-20', False, 29, [0,1], False]#90-80 NOTHING HIGH COUNTS
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-12-08', r'scan_tof_2021-12-08-17-59-46', False, 29, [0,1], False]#90-80 
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-12-08', r'scan_tof_2021-12-08-20-07-02', False, 20, [0,10], False]#90-80
    date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-07-27', r'scan_tof_2022-07-27-18-33-55', False, 50, [0,1], False]#90-80 
    
    

    


    dataPath = pl.PurePath(dataPath)
    #if not os.path.exists(os.path.join(dataPath, yearAndMonth, date, filename)):
    #    os.mkdir(os.path.join(dataPath, yearAndMonth, date, filename))
    dataFile =  os.path.join(dataPath, yearAndMonth, date, filename+r'.hdf5')

    f = FFT_ionS(filename=dataFile, dcRange=2, massRange=[5,50])
    ff=[f]

    for x in ff:
        x.read_split(overNight = True, firstTry =firstTry, sumNo = sumNo, usefulRange = usefulRange, cu=cu)
        #check the spectrum
        #for i in range(0,4):
        #   plt.plot(x.spec[i],label=str(i))
        #   np.savetxt('spec.txt',x.spec[i])
        #plt.plot(x.spec[0],label=str(0))
        #plt.legend()
        #plt.show()
        #plt.plot(np.arange(x.data.shape[0]),np.sum(x.data,1))
        plt.plot(x.calculator.pixel2mass(np.arange(x.data.shape[0])),np.sum(x.data[:,int(140/660*x.data.shape[1]):],1),label='Mass[100,600]')
        #plt.plot(np.arange(sum.shape[0]),np.sum(sum[:,int(140/660*sum.shape[1]):],1),label='Mass[100,600]')
        plt.legend()
        plt.show()
        
        x.mass_spectra()
        
        x.window(windowSize=100, direction='left')
        x.rmvExp()
        #x.smooth(windowSize=999)
        #x.rebinS(factor=20)
        #x.STFTS(windowsize=int(200/600*3600))
        x.padding(paddingSize=420000)
        #x.useFilter(1500/33*1e12,9000/33*1e12)
        x.dataProcess()
        x.show_Spectra()
        plt.show()
        x.show_FFT()
        plt.show()
        
        
        
    