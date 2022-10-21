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
    #def __init__(self, filename, massRange=[5, 70], dcRange=2, cal_mass=[43,58], cal_pixel=[345,1073]):
    def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[374,455]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[197,279]):
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
        for filename in self.filename:
            with h5py.File(filename, 'r+') as f:
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

    def read_split(self, overNight = False, firstTry = False, sumNo = 20, usefulRange = [0, 5], cu = False):
        for filename in self.filename:
            with h5py.File(filename, 'r+') as f:
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
                    self.data[filename] = np.array(f['useful'])
                    print(np.array(f['useful_description']))
                else:
                    dataSet_use = ['sum'+str(i) for i in range(usefulRange[0],usefulRange[1])]#[x for x in dataSet_created.keys()]
                    print('Now the dataset \'',  dataSet_use, '\' is in use!')
                    for d in dataSet_use:
                        self.data[filename] = self.data[filename] + np.array(f[d])
                    try:
                        f.create_dataset('useful', data=self.data[filename]) #create dataset 'sum'
                        useful_description = 'Sum of trace in range '+str((np.array(usefulRange)-[0,1])*sumNo)+' is in use!'
                        print(useful_description)
                        f.create_dataset('useful_description', data=useful_description) #create dataset 'sum'
                        print('Dataset \''+ 'useful' +'\' is created!')
                    except RuntimeError:
                        print('Dataset \''+ 'useful' +'\''+' already exists!')
                        if cu:
                            print('Begin to overwrite dataset \''+ 'useful' +'\''+'!')
                            interData = f['useful']  # load and overwrite the data
                            interData[...] = self.data[filename]

    def FFT(self, t, y):
        n = len(t)
        delta = (max(t) - min(t)) / (n-1)
        k = int(n/2)
        f = np.arange(k) / (n*delta) / 10**12  # frequency unit THz
        Y = abs(sp.fft.fft(y))[:k]
        return np.array([f, Y])

    def mass_spectra(self):
        for filename in self.filename:
            self.spectraBottle[filename]={}
            data = self.data[filename]
            for gas in self.gasBottle.keys():
                if len(self.gasBottle[gas])==2:
                    [pixelMin, pixelMax] = list(map(int, self.gasBottle[gas]))
                    self.spectraBottle[filename][gas] = np.sum(
                        data[pixelMin:pixelMax, :], 0)
                elif len(self.gasBottle[gas])==4:
                    [pixel0, pixel1, pixel2, pixel3] = list(map(int, self.gasBottle[gas]))
                    self.spectraBottle[filename][gas] = np.sum(
                        -data[pixel0:pixel1, :], 0)+np.sum(
                        data[pixel2:pixel3, :], 0)
                else:
                    [pixel0, pixel1, pixel2, pixel3, useless] = list(map(int, self.gasBottle[gas]))
                    self.spectraBottle[filename][gas] = np.sum(
                        data[pixel0:pixel1, :], 0)+np.sum(
                        data[pixel2:pixel3, :], 0)

    def window(self, windowSize, direction):
        '''
        windowSize is in fs.
        '''
        for filename in self.filename:
            for gas in self.spectraBottle[filename].keys():
                
                if direction[filename] == 'left':
                    data = self.spectraBottle[filename][gas]
                elif direction[filename] == 'right':
                    data = self.spectraBottle[filename][gas][::-1]
                if 'rebin' not in gas:
                    self.spectraBottle[filename][gas], self.delay = self.inter_window(gas, data, self.delay, windowSize=windowSize)
                else:
                    self.spectraBottle[filename][gas], self.rebin_delay = self.inter_window(gas, data, self.rebin_delay, windowSize=windowSize)
            
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
        for filename in self.filename:
            for gas in self.spectraBottle[filename].keys():
                inter_data = self.spectraBottle[filename][gas]
                if "rebin" not in gas:
                    self.spectraBottle[filename][gas], self.delay = self.inter_padding(gas, inter_data, self.delay, paddingSize = paddingSize)
                else:
                    self.spectraBottle[filename][gas], self.rebin_delay = self.inter_padding(gas, inter_data, self.rebin_delay, paddingSize = paddingSize)

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
        for filename in self.filename:
            self.fftS[filename]={}
            for gas in self.spectraBottle[filename].keys():
                if self.ifrebin:
                    if "rebin" not in gas:
                        continue
                    else:
                        delay = self.rebin_delay
                        self.fftS[filename][gas] = self.FFT(
                            delay, self.spectraBottle[filename][gas][self.windowSize::])
                else:
                    if "rebin" not in gas:
                        delay = self.delay
                        self.fftS[filename][gas] = self.FFT(
                            delay, self.spectraBottle[filename][gas][self.windowSize::])
                    else:
                        continue

    # def filterS(self, band):
     #   for gas in self.spectraBottle.keys():
      #      N, wn = sps.buttord(1e11, )
       #     fil = sps.butter(10, band, btype='band', fs=1000, output='sos')
        #    self.spectraBottle[gas] = sps.sosfilt(fil, self.spectraBottle[gas])

    def STFTS(self, legend, windowsize=0, ratio=10):
        '''
        windowsize is window size in fs
        '''
        for filename in self.filename:
            self.stftS[filename] = {}
            for gas in self.spectraBottle[filename].keys():
                #if 'rebin' not in gas or 'filter' not in gas:
                if 'rebin' not in gas:
                    continue
                f, t, self.stftS[filename][gas] = sps.stft(self.spectraBottle[filename][gas], fs=1/(
                    self.rebin_delay[2]-self.rebin_delay[1]), noverlap=windowsize-2, nperseg=windowsize, nfft=windowsize*3)
                #self.stftS[gas] = self.stftS[gas][:,::-1]*10
                self.stftS[filename][gas] = self.stftS[filename][gas]*10
                vmax=abs(self.stftS[filename][gas][:int(len(f)/15)]).max()
                vmin=abs(self.stftS[filename][gas]).min()
                norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
                levels = np.arange(vmin,vmax,(vmax-vmin)/500)
                plt.figure(figsize=(8.0, 5.0))
                plt.contourf(t*10**15+120, f/10**12*33.35641,
                            np.abs(self.stftS[filename][gas]), levels=levels, cmap='jet', norm=norm)
                #plt.xlim([100, self.delay[-1]])
                plt.ylim([2600, 4600])
                plt.ylabel('Frequency [cm-1]')
                plt.xlabel('Time [fs]')
                plt.clim(0,abs(self.stftS[filename][gas][:int(len(f)/15)]).max())
                plt.colorbar(ticks=levels[::100])
                plt.title(gas+legend[filename])
                plt.tight_layout()
                plt.savefig(str('2021-09-22+')+str(legend[filename]),dpi=720, bbox_inches='tight', pad_inches=0.2)
                plt.show()

    def rmvExp(self):
        for filename in self.filename:
            for gas in self.spectraBottle[filename].keys():
                #y = self.spectraBottle[gas]
                #k, b = np.polyfit(f.delay, np.log(y), 1)#, w=np.sqrt(y))
                #self.spectraBottle[gas] = y- np.exp(b)*np.exp(f.delay*k)
                self.spectraBottle[filename][gas]=polynomial(self.spectraBottle[filename][gas], order=10, plot=False)

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
        for filename in self.filename:
            for gas in list(self.spectraBottle[filename].keys()):
                self.spectraBottle[filename]['rebin_'+gas] = self.rebin_factor(self.spectraBottle[filename][gas], factor)
        self.rebin_delay = self.rebin_factor(self.delay, factor)

    def show_FFT(self,legend):
        fig=plt.figure()
        i=1
        for filename in self.filename:
            ax = fig.add_subplot(7,1,i)
            i=i+1
            #unused_f, ref = self.fftS['CH3OH+']
            for gas in self.fftS[filename].keys():
                [f, Y] = self.fftS[filename][gas]
                ax.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:], label=gas)
                ax.legend((gas,))
                plt.title('FFT[100:600]'+str(legend[filename]))
            ax.set_xlabel('Frequency/cm-1')
            ax.set_ylabel('a.u.')
            ax.set_ylim([0,np.max(Y[self.dcRange:1000])*3/2])
            ax.set_xlim([0,20000])
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

    def show_Spectra(self, legend, shift):
        fig = plt.figure()
        ax = fig.add_subplot()
        for filename in self.filename:
            for gas in self.spectraBottle[filename].keys():

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
                            self.spectraBottle[filename][gas]/np.sum(self.spectraBottle[filename][gas]), label=legend[filename])
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
        for filename in self.filename:
            for gas in list(self.spectraBottle[filename].keys()):
                # if 'rebin' not in gas:
                fs = 1/(self.delay[90]-self.delay[89])
                # else:
                #     fs = 1/(self.rebin_delay[90]-self.rebin_delay[89])
                self.spectraBottle[filename]['filter_'+gas] = self.butter_bandpass_filter(self.spectraBottle[filename][gas], lowcut, highcut, fs)

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
        r'C:\Users\user\Desktop\Data_newTOF\2021-09\2021-09-10')
    dataPath1 = pl.PureWindowsPath(
        r'C:\Users\user\Desktop\Data_newTOF\2021-09\2021-09-11')
    dataPath3 = pl.PureWindowsPath(
        r'C:\Users\user\Desktop\Data_newTOF\2021-09\2021-09-09')
    datafile1 =  os.path.join(dataPath1, r'scan_tof_2021-09-11-17-49-19.hdf5')#130-160
    datafile =  os.path.join(dataPath, r'scan_tof_2021-09-10-17-37-00.hdf5')#160-120
    datafile2 =  os.path.join(dataPath1, r'scan_tof_2021-09-11-14-34-38.hdf5')#110-160
    datafile3 =  os.path.join(dataPath1, r'scan_tof_2021-09-11-13-20-25.hdf5')#90-160
    datafile4 =  os.path.join(dataPath1, r'scan_tof_2021-09-11-11-45-33.hdf5') #70-160
    datafile5 =  os.path.join(dataPath1, r'scan_tof_2021-09-11-10-28-14.hdf5')#60-160
    datafile6 =  os.path.join(dataPath1, r'scan_tof_2021-09-11-19-12-43.hdf5')#120-180

    datafile_long_1 = os.path.join(dataPath3, r'scan_tof_2021-09-09-21-37-12.hdf5')
    
    
    legend = {datafile1:'130-160',datafile:'160-120', datafile2:'110-160', datafile3:'90-160', datafile4:'70-160', datafile5:'60-120', datafile6:'120-180',\
        datafile_long_1:'90-120'}
    shift = {datafile1:0,datafile:0,datafile2:0,datafile3:0,datafile4:0,datafile5:0,datafile6:0,datafile_long_1:0}

    #plt.plot(delay*2*3.33564,np.sum(sum[int(cal.cal_pixel(30.5)):int(cal.cal_pixel(31.5))],0),label='Mass 32')
    
    #check the spectrum
    #for i in range(0,20):
    #    plt.plot(spec[i],label=str(i))
    #plt.plot(spec[0],label=str(0))
    #plt.legend()
    #plt.show()
    #print('Intensity changes'+str(np.sum(np.mean(spec[20:40],0))/np.sum(spec[0])))


    #f = FFT_ionS(filename=[datafile1,datafile2,datafile3,datafile5], dcRange=2, massRange=[5,50])
    f = FFT_ionS(filename=[datafile_long_1], dcRange=2, massRange=[5,50])
    ff=[f]

    for x in ff:
        x.read_split(overNight = True, firstTry = False, sumNo = 20, usefulRange = [0,1], cu=False)
        xdirection = {datafile1:'left', datafile:'right', datafile2:'left', datafile3:'left', datafile4:'left', datafile5:'left', datafile6:'left',datafile_long_1:'right'}
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

        x.mass_spectra()
        
        x.window(windowSize=150, direction=xdirection)
        x.rmvExp()
        #x.smooth(windowSize=301)
        #x.rebinS(factor=20)
        #x.STFTS(legend=legend,windowsize=int(200/600*3500))
        #x.smooth(windowSize=199)
        x.padding(paddingSize=140000)
        #x.useFilter(100/33*1e12,5000/33*1e12) 
        x.dataProcess()
        x.show_FFT(legend=legend)
        x.show_Spectra(legend=legend,shift=shift)
        plt.show()
        
        
        
    