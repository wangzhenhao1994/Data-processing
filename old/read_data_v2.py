import sys
import os
import pathlib as pl
from matplotlib.colors import Normalize
from matplotlib import cm
from calculate_k_b import Calibration_mass
import numpy as np
import scipy as sp
import scipy.signal as sps
import scipy.interpolate as spi
from obspy.signal.detrend import polynomial
import matplotlib as mpl
from matplotlib import pyplot as plt
import h5py

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
    def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[433,514]):

        '''

        '''
        self.filename = filename
        self.channelSize = 1024
        self.scanLength = 70000#72000
        self.peakRange = [-15, 15]  # range of the peak
        self.delay = np.arange(self.scanLength)/self.scanLength*(116.7/120)*94*2*3.33564*10**-15
        self.delay = self.delay - 60.00*10**-15
        self.rebin_delay = None
        self.ifrebin = False
        self.longstage = True
        self.data = {}
        self.massRange = massRange
        self.dcRange = dcRange
        self.windowSize = 100
        self.calculator = Calibration_mass(cal_mass, cal_pixel)
        self.gasBottle = {
            #"43": self.calculator.cal_pixel(43) + self.peakRange,
            #"45": self.calculator.cal_pixel(45) + self.peakRange,
            #"58": self.calculator.cal_pixel(58) + self.peakRange,
            #"CH3OH+": self.calculator.cal_pixel(32) + self.peakRange,
            #"CH2OH+": self.calculator.cal_pixel(31) + self.peakRange,
            #"CHOH+": self.calculator.cal_pixel(30) + self.peakRange,
            #"CHO+": self.calculator.cal_pixel(29) + self.peakRange,
            #"CO+":self.calculator.cal_pixel(28) +self.peakRange,
            #"CH2+":self.calculator.cal_pixel(14) +self.peakRange,
            #"CH3+":self.calculator.cal_pixel(15) +self.peakRange,
            # "N": self.calculator.cal_pixel(14) +self.peakRange,
            # "O": self.calculator.cal_pixel(16)+self.peakRange,
            #"H2O+": self.calculator.cal_pixel(18)+self.peakRange,
            # "N2": self.calculator.cal_pixel(28)+self.peakRange,
            # "$O_2$": self.calculator.cal_pixel(32)+self.peakRange,
            # "N": self.calculator.cal_pixel(14) +self .peakRange,
            # "O": self.calculator.cal_pixel(16)+self.peakRange,
            "H2O+": self.calculator.cal_pixel(18)+self.peakRange,
            #"OH-": self.calculator.cal_pixel(17)+self.peakRange,
            # "N2": self.calculator.cal_pixel(28)+self.peakRange,
            # "O2": self.calculator.cal_pixel(32)+self.peakRange,
            # "Ar": self.calculator.cal_pixel(32)+self.peakRange,
            #"\'CH3OH+\'-\'CH2OH+\'":np.append(self.calculator.cal_pixel(43) + self.peakRange, self.calculator.cal_pixel(58) + self.peakRange),
            #"\'CH3OH+\'+\'CH2OH+\'": #[0,self.channelSize]
            #np.append(np.append(self.calculator.cal_pixel(17) + self.peakRange, self.calculator.cal_pixel(18) + self.peakRange),0),
        }
        self.spectraBottle = {}
        self.fftS = {}
        self.stftS = {}
        self.spec = None
        self.scanNo = None

    def checkData(self):
        for key in self.filename.keys():
            with h5py.File(self.filename[key], 'r+') as f:
                print(np.array(f['parameters']))
                bigger = True
                i = 0
                while bigger:
                    try:
                        a = f['data'][i*self.channelSize,0]
                        i+=1
                    except ValueError:
                        bigger = False
                        self.scanNo = i
                        f.create_dataset('scanNo', data=i)
                        print('The number of scan is '+str(i)+'!')

    def checkPrPuCount(self):
        for key in self.filename.keys():
            with h5py.File(self.filename[key], 'r+') as f:
                
                ppNum = np.array(f['ppNum'])
                #print(ppNum)
                prOnly = 0
                puOnly = 0
                pupr = 0
                prOnly=np.sum(np.array(f['data'][int(10.3*self.channelSize):(10+1)*self.channelSize]))
                
                puOnly=np.sum(np.array(f['data'][int(11.3*self.channelSize):(11+1)*self.channelSize]))
                pupr=np.sum(np.array(f['data'][int(12.3*self.channelSize):(12+1)*self.channelSize]))
                #plt.plot(puOnly)
                #plt.show()
                diff = (prOnly+puOnly-pupr)/pupr*100
                print(diff)

    def read_split(self, overNight = False, firstTry = False, sumNo = 20, usefulRange = [0, 5], cu = False):
        for key in self.filename.keys():
            with h5py.File(self.filename[key], 'r+') as f:
                print(f.keys())
                self.spec=np.array(f['spectrumLog'])
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
                                    np.array(f['data'][i*self.channelSize:(i+1)*self.channelSize])
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
                    self.data[key] = np.array(f['useful'])
                    print(np.array(f['useful_description']))
                else:
                    dataSet_use = ['sum'+str(i) for i in range(usefulRange[0],usefulRange[1])]#[x for x in dataSet_created.keys()]
                    print('Now the dataset \'',  dataSet_use, '\' is in use!')
                    for d in dataSet_use:
                        self.data[key] = self.data[key] + np.array(f[d])
                    try:
                        f.create_dataset('useful', data=self.data[key]) #create dataset 'sum'
                        useful_description = 'Sum of trace in range '+str((np.array(usefulRange)-[0,1])*sumNo)+' is in use!'
                        print(useful_description)
                        f.create_dataset('useful_description', data=useful_description) #create dataset 'sum'
                        print('Dataset \''+ 'useful' +'\' is created!')
                    except RuntimeError:
                        print('Dataset \''+ 'useful' +'\''+' already exists!')
                        if cu:
                            print('Begin to overwrite dataset \''+ 'useful' +'\''+'!')
                            interData = f['useful']  # load and overwrite the data
                            interData[...] = self.data[key]

    def FFT(self, t, y):
        n = len(t)
        delta = (max(t) - min(t)) / (n-1)
        k = int(n/2)
        f = np.arange(k) / (n*delta) / 10**12  # frequency unit THz
        Y = abs(sp.fft.fft(y))[:k]
        return np.array([f, Y])

    def mass_spectra(self):
        for key in self.filename.keys():
            self.spectraBottle[key]={}
            data = self.data[key]
            for gas in self.gasBottle.keys():
                if len(self.gasBottle[gas])==2:
                    [pixelMin, pixelMax] = list(map(int, self.gasBottle[gas]))
                    self.spectraBottle[key][gas] = np.sum(
                        data[pixelMin:pixelMax, :], 0)
                elif len(self.gasBottle[gas])==4:
                    [pixel0, pixel1, pixel2, pixel3] = list(map(int, self.gasBottle[gas]))
                    self.spectraBottle[key][gas] = np.sum(
                        -data[pixel0:pixel1, :], 0)+np.sum(
                        data[pixel2:pixel3, :], 0)
                else:
                    [pixel0, pixel1, pixel2, pixel3, useless] = list(map(int, self.gasBottle[gas]))
                    self.spectraBottle[key][gas] = np.sum(
                        data[pixel0:pixel1, :], 0)+np.sum(
                        data[pixel2:pixel3, :], 0)

    def window(self, windowSize, direction):
        '''
        windowSize is in fs.
        '''
        for key in self.filename.keys():
            for gas in self.spectraBottle[key].keys():
                
                if direction[filename] == 'left':
                    data = self.spectraBottle[key][gas]
                elif direction[filename] == 'right':
                    data = self.spectraBottle[key][gas][::-1]
                if 'rebin' not in gas:
                    self.spectraBottle[key][gas], self.delay = self.inter_window(gas, data, self.delay, windowSize=windowSize)
                else:
                    self.spectraBottle[key][gas], self.rebin_delay = self.inter_window(gas, data, self.rebin_delay, windowSize=windowSize)
            
    def inter_window(self, gas, data, delay, windowSize=0, direction='left'):
        '''
        windowSize is in fs.
        '''
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
        for key in self.filename.keys():
            for gas in self.spectraBottle[key].keys():
                inter_data = self.spectraBottle[key][gas]
                if "rebin" not in gas:
                    self.spectraBottle[key][gas], self.delay = self.inter_padding(gas, inter_data, self.delay, paddingSize = paddingSize)
                else:
                    self.spectraBottle[key][gas], self.rebin_delay = self.inter_padding(gas, inter_data, self.rebin_delay, paddingSize = paddingSize)

    def inter_padding(self, gas, inter_data, delay, paddingSize=0):
        data = None
        delayStep = delay[1]-delay[0]
        delay = np.concatenate((
            np.arange(delay[0]-(paddingSize+1) *
                      delayStep, delay[0], delayStep),
            delay,
            np.arange(delay[-1]+delayStep, delay[-1] +
                      (paddingSize)*delayStep, delayStep)
        ))
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
        for key in self.filename.keys():
            self.fftS[key]={}
            for gas in self.spectraBottle[key].keys():
                if self.ifrebin:
                    if "rebin" not in gas:
                        continue
                    else:
                        delay = self.rebin_delay
                        self.fftS[key][gas] = self.FFT(
                            delay, self.spectraBottle[key][gas][self.windowSize::])
                else:
                    if "rebin" not in gas:
                        delay = self.delay
                        self.fftS[key][gas] = self.FFT(
                            delay, self.spectraBottle[key][gas][self.windowSize::])
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
        for key in self.filename.keys():
            self.stftS[key] = {}
            for gas in self.spectraBottle[key].keys():
                #if 'rebin' not in gas or 'filter' not in gas:
                if 'rebin' not in gas:
                    continue
                f, t, self.stftS[key][gas] = sps.stft(self.spectraBottle[key][gas], fs=1/(
                    self.rebin_delay[2]-self.rebin_delay[1]), noverlap=windowsize-2, nperseg=windowsize, nfft=windowsize*3)
                #self.stftS[gas] = self.stftS[gas][:,::-1]*10
                self.stftS[key][gas] = self.stftS[key][gas]*10
                vmax=abs(self.stftS[key][gas][:int(len(f)/15)]).max()
                vmin=abs(self.stftS[key][gas]).min()
                norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
                levels = np.arange(vmin,vmax,(vmax-vmin)/500)
                plt.figure(figsize=(8.0, 5.0))
                plt.contourf(t*10**15+120, f/10**12*33.35641,
                            np.abs(self.stftS[key][gas]), levels=levels, cmap='jet', norm=norm)
                #plt.xlim([100, self.delay[-1]])
                plt.ylim([2600, 4600])
                plt.ylabel('Frequency [cm-1]')
                plt.xlabel('Time [fs]')
                plt.clim(0,abs(self.stftS[key][gas][:int(len(f)/15)]).max())
                plt.colorbar(ticks=levels[::100])
                plt.title(gas+key)
                plt.tight_layout()
                plt.savefig('2021-10-18+'+key,dpi=720, bbox_inches='tight', pad_inches=0.2)
                #plt.show()

    def rmvExp(self):
        for key in self.filename.keys():
            for gas in self.spectraBottle[key].keys():
                #y = self.spectraBottle[gas]
                #k, b = np.polyfit(f.delay, np.log(y), 1)#, w=np.sqrt(y))
                #self.spectraBottle[gas] = y- np.exp(b)*np.exp(f.delay*k)
                self.spectraBottle[key][gas]=polynomial(self.spectraBottle[key][gas], order=10, plot=False)

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
        for key in self.filename.keys():
            for gas in list(self.spectraBottle[key].keys()):
                self.spectraBottle[key]['rebin_'+gas] = self.rebin_factor(self.spectraBottle[key][gas], factor)
        self.rebin_delay = self.rebin_factor(self.delay, factor)

    def show_FFT(self):
        fig, axx = plt.subplots(3, 2, sharey=True, sharex=True)
        ax = [axx[i,j] for i in range(3) for j in range(2)]
        print(ax.__len__())
        i=0
        for key in self.filename.keys():
            #unused_f, ref = self.fftS['CH3OH+']
            for gas in self.fftS[key].keys():
                [f, Y] = self.fftS[key][gas]
                ax[i].plot(f[self.dcRange:]*33.35641, Y[self.dcRange:], label=gas+':'+key)
                ax[i].legend()
                #ax[i].set_title('FFT[100:600]'+key)
                ax[i].set_xlabel('Frequency/cm-1')
                ax[i].set_ylabel('a.u.')
                ax[i].set_ylim([0,np.max(Y[self.dcRange:1000])*3/2])
                ax[i].set_xlim([0,5000])
            i=i+1
            #plt.tight_layout()
            #plt.savefig("fft_31"+d,dpi=720, bbox_inches='tight', pad_inches=0.2)

            # plt.figure()
            # for gas in ['N', 'O']:
            #    [f, Y] = self.fftS[gas]
            #    plt.plot(f[self.dcRange:], Y[self.dcRange:], label=gas)
            #    plt.xlabel('THz')
            #    plt.ylabel('a.b.')
            #    plt.legend(gas)
            # plt.legend()
        plt.show()

    def show_Spectra(self, shift):
        fig = plt.figure()
        ax = fig.add_subplot()
        for key in self.filename.keys():
            for gas in self.spectraBottle[key].keys():
                #ax = fig.add_subplot(2,1,i)
                #if 'rebin' not in gas: #or 'filter' not in gas:
                if 'filter' not in gas:
                   continue
                #elif 'rebin' not in gas:
                    #delay = self.delay
                else:
                    #delay = self.rebin_delay
                    delay = self.delay
                    ax.plot(delay*10**15+shift[filename],
                            self.spectraBottle[key][gas]/np.sum(self.spectraBottle[key][gas]), label=key)
                #plt.plot(self.spectraBottle[gas]/np.sum(self.spectraBottle[gas]), label=gas)
        ax.set_xlabel("Delay/fs")
        ax.set_ylabel('a.u.')
        ax.legend()
        #plt.xlim([200,600])
        
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
        for key in self.filename.keys():
            for gas in list(self.spectraBottle[key].keys()):
                # if 'rebin' not in gas:
                fs = 1/(self.delay[90]-self.delay[89])
                # else:
                #     fs = 1/(self.rebin_delay[90]-self.rebin_delay[89])
                self.spectraBottle[key]['filter_'+gas] = self.butter_bandpass_filter(self.spectraBottle[key][gas], lowcut, highcut, fs)

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
        dataPath = r'C:\Users\user\Desktop\Data_newTOF'
        dataPath = pl.PureWindowsPath(dataPath)
    yearAndMonth = r'2021-10'
    
    specDict = {
        '120-120': [r'2021-10-16', r'scan_tof_2021-10-16-19-47-59', False, 85, [0,1], False, 'left'],
        '100-110': [r'2021-10-14', r'scan_tof_2021-10-14-21-00-11', False, 36, [0,1], False, 'left'],
        #'90-110' :[r'2021-10-16', r'scan_tof_2021-10-16-16-57-38', False, 36, [0,1], False, 'left'],#90-110 35 Ar, no signal
        '90-110':  [r'2021-10-14', r'scan_tof_2021-10-14-19-17-24', False, 36, [0,1], False, 'left'],
        '80-110':  [r'2021-10-14', r'scan_tof_2021-10-14-17-50-39', False, 36, [0,1], False, 'left'],
        '70-110':  [r'2021-10-14', r'scan_tof_2021-10-14-16-08-58', False, 36, [0,1], False, 'left'],
        '60-110':  [r'2021-10-14', r'scan_tof_2021-10-14-14-40-33', False, 36, [0,1], False, 'left'],

        '120-90' : [r'2021-10-15', r'scan_tof_2021-10-15-11-22-53', False, 36, [0,1], False, 'left'],#120-90 35 Ar # two CEP, try again #no signal
        '110-90' : [r'2021-10-15', r'scan_tof_2021-10-15-14-57-49', False, 36, [0,1], False, 'left'],#110-90 35 Ar# two CEP, try again #no signal
        #'110-90' :[r'2021-10-15', r'scan_tof_2021-10-15-13-00-34', False, 36, [0,1], False, 'left'],#110-90 35 Ar# two CEP, try again #no signal
        '100-90' : [r'2021-10-14', r'scan_tof_2021-10-14-13-05-18', False, 36, [0,1], False, 'left'],#100-90 35 Ar #no signal
        '90-90':  [r'2021-10-14', r'scan_tof_2021-10-14-11-20-21', False, 36, [0,1], False, 'left'],#90-90 35 Ar #no signal
        '85-90':  [r'2021-10-15', r'scan_tof_2021-10-15-18-41-16', False, 36, [0,1], False, 'left'],#85-90 35 Ar #no signal
        '80-90':  [r'2021-10-14', r'scan_tof_2021-10-14-09-48-02', False, 36, [0,1], False, 'left'],#80-90 35 Ar #spectra has blue drift,try again later #no signal
        '75-90':  [r'2021-10-15', r'scan_tof_2021-10-15-16-56-57', False, 36, [0,1], False, 'left'],#75-90 35 Ar, two peaks
        '70-90':  [r'2021-10-13', r'scan_tof_2021-10-13-19-52-07', False, 36, [0,1], False, 'left'],#70-90 35 Ar
        '60-90':  [r'2021-10-13', r'scan_tof_2021-10-13-18-12-42', False, 36, [0,1], False, 'left'],#60-90 35 Ar

        '100-100': [r'2021-10-13', r'scan_tof_2021-10-13-16-36-09', False, 36, [0,1], False, 'left'],#100-100 35 Ar try again later, saturation.
        '90-100' : [r'2021-10-12', r'scan_tof_2021-10-13-12-49-40', False, 36, [0,1], False, 'left'],#90-100 35 Ar
        '85-100' : [r'2021-10-15', r'scan_tof_2021-10-15-20-58-30', False, 36, [0,1], False, 'left'],#85-100 35 Ar, no signal
        #'80-100'  :[r'2021-10-13', r'scan_tof_2021-10-13-14-53-44', False, 36, [0,1], Fals, 'left'e],#80-100 35 Ar, strong signal
        '80-100' : [r'2021-10-12', r'scan_tof_2021-10-12-20-17-16', False, 36, [0,1], False, 'left'],#80-100 35 Ar, strong signal
        #'80-100'  :[r'2021-10-16', r'scan_tof_2021-10-16-13-15-03', False, 36, [0,1], Fals, 'left'e],#80-100 35 Ar, no signal
        #'80-100'  :[r'2021-10-16', r'scan_tof_2021-10-16-15-19-13', True, 36, [0,1], True], 'left',  #80-100 35 Ar, little signal
        '75-100' : [r'2021-10-16', r'scan_tof_2021-10-16-11-29-18', False, 36, [0,1], False, 'left'],#75-100 35 Ar
        '70-100' : [r'2021-10-13', r'scan_tof_2021-10-13-11-05-44', False, 36, [0,1], False, 'left'],#70-100 35 Ar
        #'60-100'  :[r'2021-10-13', r'scan_tof_2021-10-13-09-20-24', False, 36, [0,1], Fals, 'left'e],#60-100 35 Ar
        '60-100' : [r'2021-10-12', r'scan_tof_2021-10-12-21-58-45', False, 15, [0,4], False, 'left'],#60-100 250 Ar
        '60-80'  : [r'2021-10-13', r'scan_tof_2021-10-13-21-40-53', False, 15, [0,5], False, 'left'],#60-80 250 Ar #no signal
        '110-70' : [r'2021-10-15', r'scan_tof_2021-10-15-09-39-43', False, 36, [0,1], False, 'left'],#110-70 35 Ar #no signal # two CEP, try again
        '110-60' : [r'2021-10-14', r'scan_tof_2021-10-14-22-57-02', False, 15, [0,4], False, 'left']#110-60 150 Ar #no signal
    }  

    pump_intensity = [60,70,80,90,100]
    probe_intensity = [90]#[90,100,110]

    xdirection = {}
    dataFile = {}
    shift = {}
    dataPath = pl.PurePath(dataPath)
    for prIntensity in probe_intensity:
        for puIntensity in pump_intensity:
            leg = str(puIntensity)+'-'+str(prIntensity)
            s = 0
            date, filename, firstTry, sumNo, usefulRange, cu, direction = specDict[leg]
            xdirection[leg]=direction
            shift[leg] = s
            file =  os.path.join(dataPath, yearAndMonth, date, filename+r'.hdf5')
            dataFile[leg] = file

    #if not os.path.exists(os.path.join(dataPath, yearAndMonth, date, filename)):
    #    os.mkdir(os.path.join(dataPath, yearAndMonth, date, filename))

    f = FFT_ionS(filename=dataFile, dcRange=2, massRange=[5,50])
    ff=[f]


    for x in ff:
        #x.read_split(overNight = True, firstTry =False, sumNo = 35, usefulRange = [0,1], cu=False)
        #check the spectrum
        #for i in range(0,5):
        #   plt.plot(x.spec[i],label=str(i))
        #plt.plot(x.spec[0],label=str(0))
        #plt.legend()
        #plt.show()
        for filename in x.filename:
            #plt.plot(np.arange(x.data[filename].shape[0]),np.sum(x.data[filename],1))
        #plt.plot(x.calculator.pixel2mass(np.arange(x.data.shape[0])),np.sum(x.data[:,int(140/660*x.data.shape[1]):],1))
        #plt.plot(np.arange(sum.shape[0]),np.sum(sum[:,int(140/660*sum.shape[1]):],1),label='Mass[100,600]')
        #plt.legend()
            plt.show()

        #x.mass_spectra()
        
        #x.window(windowSize=120, direction=xdirection)
        x.checkPrPuCount()
        #x.rmvExp()
        #x.smooth(windowSize=301)
        #x.rebinS(factor=20)
        #x.STFTS(windowsize=int(200/600*3500))
        #x.smooth(windowSize=199)
        #x.padding(paddingSize=140000)
        #x.useFilter(100/33*1e12,5000/33*1e12) 
        #x.dataProcess()
        #x.show_FFT()
        #x.show_Spectra(shift=shift)
        #plt.show()
        
        
        
    