from operator import xor
from tkinter.messagebox import YES

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
mpl.style.use('seaborn-poster')
import h5py
import pathlib as pl
import os
from decimal import Decimal
#from pynufft import NUFFT

my_path = os.path.abspath(__file__)


mpl.rcParams['lines.linewidth'] = 1
plt.rcParams["font.family"] = "arial"
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.it'] = 'Arial:italic'
mpl.rcParams['mathtext.rm'] = 'Arial'

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

    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[60,61]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[1,2], cal_pixel=[100,244]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[17,18], cal_pixel=[197,279]):
    #def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[2,16], cal_pixel=[387,1281]):
    def __init__(self, filename, massRange=[5, 50], dcRange=2, cal_mass=[2,16], cal_pixel=[3087,10241]):
        '''

        '''
        
        self.filename = filename
        self.channelSize = 12032#24000#24000#1536
        self.scanLength = 3320#1200
        self.peakRange = [-100, 100]  # range of the peak
        self.delay = np.arange(self.scanLength)/self.scanLength*(11/11.2)*100*2*2*3.33564*10**-15
        self.rebin_delay = None
        self.ifrebin = False
        self.longstage = True
        self.data = 0
        self.massRange = massRange
        self.dcRange = dcRange
        self.windowSize = 100
        self.calculator = Calibration_mass(cal_mass, cal_pixel)
        self.gasBottle = {
            #"8": self.calculator.cal_pixel(8) + self.peakRange,
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
            #"O+": self.calculator.cal_pixel(16)+self.peakRange,
            #"H2O+": self.calculator.cal_pixel(18)+self.peakRange,
            # "N2": self.calculator.cal_pixel(28)+self.peakRange,
            # "$O_2$": self.calculator.cal_pixel(32)+self.peakRange,
            # "N": self.calculator.cal_pixel(14) +self .peakRange,
            "O+": self.calculator.cal_pixel(16)+self.peakRange,
            #"H2O+": self.calculator.cal_pixel(18)+self.peakRange,
            #"OH+": self.calculator.cal_pixel(17)+self.peakRange,
            # "N2": self.calculator.cal_pixel(28)+self.peakRange,
            # "O2": self.calculator.cal_pixel(32)+self.peakRange,
            # "Ar": self.calculator.cal_pixel(32)+self.peakRange,
            "H+": self.calculator.cal_pixel(1)+self.peakRange,
            "H2+": self.calculator.cal_pixel(2)+self.peakRange,
            #"O+H2": #[0,self.channelSize]
            #np.append(self.calculator.cal_pixel(16) + self.peakRange, self.calculator.cal_pixel(2) + self.peakRange),
            #"\'15+\'  +  \'H+\'": #[0,self.channelSize]
            #np.append(self.calculator.cal_pixel(16) + self.peakRange, self.calculator.cal_pixel(15) + self.peakRange),
            #"O-H2": #[0,self.channelSize]
            #np.append(np.append(self.calculator.cal_pixel(16) + self.peakRange, self.calculator.cal_pixel(2) + self.peakRange), 0),
            #'test':self.calculator.cal_pixel(1.5)+self.peakRange,
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
        self.molecule = 'H2O'
        self.intensity = 5
        self.stage = 'piezo'
        self.exportPath = pl.PurePath(r'C:\Users\user\Desktop\Data_newTOF\plotdata\20220516')

    def checkData(self):
        for filename in self.filename:
            with h5py.File(filename, 'r+') as f:
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
    def read_multi(self):
        for file in self.filename:
            with h5py.File(file, 'r+') as f:
                self.data=self.data+np.array(f['useful'])

    def read_split(self, overNight = False, firstTry = False, sumNo = 20, usefulRange = [0, 5], cu = False):
        for filename in self.filename:
            with h5py.File(filename, 'r+') as f:
                print(np.array(f['dataD'].shape))
                plt.plot(np.sum(np.array(f['dataD'][0*self.channelSize:(1)*self.channelSize]),0))
                #print('the pump-only count is '+str(np.sum(np.array(f['dataD'][164*self.channelSize:(164+1)*self.channelSize]))/116.7),'\nthe probe-only count is '+str(np.sum(np.array(f['dataD'][165*self.channelSize:(165+1)*self.channelSize]))/116.7) )
                plt.show()
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
                    dataSet_created = {'sum0':[0,self.scanNo-1]}
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
                #print(pixelMin, pixelMax)
                inter= np.sum(
                    self.data[pixelMin:pixelMax, :], 0)
                #print(pixelMin, pixelMax)
                self.spectraBottle[gas] = inter
            elif len(self.gasBottle[gas])==4:
                [pixel0, pixel1, pixel2, pixel3] = list(map(int, self.gasBottle[gas]))
                inter1=np.sum(self.data[pixel0:pixel1, :], 0)
                inter2=np.sum(self.data[pixel2:pixel3, :], 0)
                self.spectraBottle[gas] =inter1 + inter2
                print('Sum two mass!!')
            elif len(self.gasBottle[gas])==5:
                [pixel0, pixel1, pixel2, pixel3, useless] = list(map(int, self.gasBottle[gas]))
                inter1=np.sum(self.data[pixel0:pixel1, :], 0)
                inter2=np.sum(self.data[pixel2:pixel3, :], 0)
                self.spectraBottle[gas] =inter1 - inter2
                print('Sub two mass!!')


    def window(self, windowSize=0, direction='left', shift=0):
        '''
        windowSize is in fs.
        '''
        for gas in self.spectraBottle.keys():
            if 'rebin' not in gas:
                self.spectraBottle[gas], self.delay = self.inter_window(gas, self.delay, windowSize=windowSize, direction=direction, shift=shift)
            else:
                self.spectraBottle[gas], self.rebin_delay = self.inter_window(gas, self.rebin_delay, windowSize=windowSize, direction=direction)
            
    def inter_window(self, gas, delay, windowSize=0, direction='left',shift=0):
        '''
        windowSize is in fs.
        '''
        data = self.spectraBottle[gas]
        __len = np.size(data)
        windowSize = int(windowSize/(delay[1]-delay[0])*1e-15) 
        if direction == 'left':
            if gas == "H2+":
                n = int(shift/(delay[1]-delay[0])*1e-15)
                window = data[-(__len-windowSize+n+1):-n-1]
                #print(n)
            else:
                window = data[-(__len-windowSize+1):-1]
            #window = data[-(__len-windowSize):]
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
        window=polynomial(window, order=3, plot=False)
        hwindow = np.hanning(len(window))#sps.windows.flattop(len(window),sym=FALSE)
        window2=window#*hwindow
        return window2, delay

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
            #np.arange(delay[0]-(paddingSize+1) *
            #          delayStep, delay[0], delayStep),
            delay,
            np.arange(delay[-1]+delayStep, delay[-1] +
                      (paddingSize)*delayStep, delayStep)
        ))
        inter_data = self.spectraBottle[gas]
        data = np.concatenate(
            #(np.zeros(paddingSize)+inter_data[0], inter_data, np.zeros(paddingSize)+inter_data[-1]), axis=0)
            (inter_data, (np.zeros(paddingSize))), axis=0)
        delay = delay[:len(data)]
        #data=data[:len(delay)]
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
            vmax=abs(self.stftS[gas][:int(len(f)/15)]).max()
            vmin=abs(self.stftS[gas]).min()
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            levels = np.arange(vmin,vmax,(vmax-vmin)/500)
            plt.figure(figsize=(8.0, 5.0))
            plt.contourf(t*10**15+120, f/10**12*33.35641,
                           np.abs(self.stftS[gas]), levels=levels, cmap='jet', norm=norm)
            #plt.xlim([100, self.delay[-1]])
            plt.ylim([0, 5000])
            plt.ylabel('Frequency [cm-1]')
            plt.xlabel('Time [fs]')
            plt.clim(0,abs(self.stftS[gas][:int(len(f)/15)]).max())
            plt.colorbar(ticks=levels[::100])
            d = '_100_120'
            plt.title(gas+d)
            plt.tight_layout()
            plt.savefig("20211012_sft_18"+d,dpi=720, bbox_inches='tight', pad_inches=0.2)
            plt.show()

    def rmvExp(self):
        for gas in self.spectraBottle.keys():
            #y = self.spectraBottle[gas]
            #k, b = np.polyfit(f.delay, np.log(y), 1)#, w=np.sqrt(y))
            #self.spectraBottle[gas] = y- np.exp(b)*np.exp(f.delay*k)
            self.spectraBottle[gas]=polynomial(self.spectraBottle[gas], order=3, plot=False)

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
            #if 'filter' not in gas:
            #    continue
            [f, Y] = self.fftS[gas]
            #Y=np.where(np.logical_and(f>=10, f<=4500),Y,0)
            if gas == 'H2+':
                plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=r'$\mathrm{H_2^+}$')
            if gas == 'O+':
                plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=r'$\mathrm{O^+}$')
            if gas == 'O+H2':
                plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=r'$\mathrm{H_2^++O^+}$')
                print("The max of the FFT is ",np.amax( Y[self.dcRange:]/1000))
            if gas == 'O-H2':
                plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=r'$\mathrm{H_2^+-O^+}$')
            else:
                plt.plot(f[self.dcRange:]*33.35641, Y[self.dcRange:]/1000, label=gas)
                pass
            plt.xlabel(r'$\mathrm{Frequency (cm^{-1})}$')
            plt.ylabel('a.u.')
            plt.legend(loc=2)
            #plt.xlabel('Frequency/cm-1')
            #plt.ylim([0,np.max(Y[self.dcRange:1000])*3/2])
            plt.xlim([0,4500])
            np.savetxt(pl.PureWindowsPath(self.exportPath, str('FFT')+str(gas)+self.molecule+str('%.1E' % Decimal(self.intensity))+str('.dat')), np.abs(Y))
            np.savetxt(pl.PureWindowsPath(self.exportPath, str('FFT')+str(gas)+self.molecule+str('%.1E' % Decimal(self.intensity))+str('frequency')+str('.dat')), f*33.35641)
        #plt.legend()
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
        #plt.close()
        return np.amax( Y[self.dcRange:]/1000)

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
            if gas == 'H2+':
                plt.plot(delay*10**15+shift,self.spectraBottle[gas],label=r'$\mathrm{H_2^+}$')
            if gas == 'O+':
                plt.plot(delay*10**15+shift,self.spectraBottle[gas],label=r'$\mathrm{O^+}$')
            if gas == 'O+H2':
                plt.plot(delay*10**15+shift,self.spectraBottle[gas],label=r'$\mathrm{H_2^++O^+}$')
            if gas == 'O-H2':
                plt.plot(delay*10**15+shift,self.spectraBottle[gas],label=r'$\mathrm{H_2^+-O^+}$')
            
            #plt.plot(self.spectraBottle[gas]/np.sum(self.spectraBottle[gas]), label=gas)
            plt.xlabel("Delay (fs)")
            plt.ylabel('Count (a.u.)')
            #plt.xlim([200,600])
            plt.legend()
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
            fs = 1/(self.delay[1]-self.delay[0])
            # else:
            #     fs = 1/(self.rebin_delay[90]-self.rebin_delay[89])
            self.spectraBottle['filter_'+gas] = self.butter_bandpass_filter(self.spectraBottle[gas], lowcut, highcut, fs)

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

class FFTOSC(FFT_ionS):
    def __init__(self, delay, data, massRange=[31, 32], dcRange=2, label=None):
        super().__init__(delay, data, massRange=massRange, dcRange=dcRange, label=label)

    def mass_spectra(self):
        self.spectraBottle['CH3OH+'] = self.data


if __name__ == '__main__':
    
    dataPath = r'D:\DataProcessing'
    yearAndMonth = r'2022-10'
    
    
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-12', r'scan_tof_2021-10-12-20-17-16', False, 36, [0,1], False]#80-100 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-13', r'scan_tof_2021-10-13-11-05-44', False, 36, [0,1], False]#70-100 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-13', r'scan_tof_2021-10-13-09-20-24', False, 36, [0,1], False]#60-100 35 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-12', r'scan_tof_2021-10-12-21-58-45', False, 15, [0,4], True]#60-100 250 Ar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2021-10-12', r'scan_tof_2021-10-12-16-05-08', False, 30, [0,1], False]

    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-08', r'scan_tof_2022-04-08-00-45-34', False, 10, [0,1], True] #200-200
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-08', r'scan_tof_2022-04-08-19-51-32', False, 10, [0,1], True] #160-160
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-09', r'scan_tof_2022-04-09-17-34-37', False, 10, [0,1], True] #120-120

    #2022-04-26
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-26', r'scan_tof_2022-04-26-16-47-41', False, 10, [0,1], True] #160-160 MCP 3700V
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-26', r'scan_tof_2022-04-26-18-38-32', False, 10, [0,1], True] #160-160 MCP 3400V #H+ and H2+ and O+ count too little, oscillation in H2O+
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-27', r'scan_tof_2022-04-27-10-54-39', False, 10, [0,1], True] #160-160 MCP 2900V
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-27', r'scan_tof_2022-04-27-12-10-32', False, 10, [0,1], True] #160-160 MCP 3200V
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-27', r'scan_tof_2022-04-27-14-15-19', False, 10, [0,1], True] #160-160 MCP 3200V
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-28', r'scan_tof_2022-04-28-12-07-57', False, 10, [0,1], True] #160-160 MCP 3400V 8e-8mbar
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-28', r'scan_tof_2022-04-28-19-19-47 - Copy', False, 190, [0,1], True] #160-160 MCP 3000V 1.5e-8mbar 1 bind width, only MAss 18
    
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-29', [r'scan_tof_2022-04-29-11-33-27'], False, 190, [0,1], True] #160-160 MCP 3220V Mass 1 and 17
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-29', [r'scan_tof_2022-04-29-14-23-07'], False, 190, [0,1], True] #160-160 MCP 3220V Mass 1 and 17
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-29', [r'scan_tof_2022-04-29-16-18-15'], False, 190, [0,1], True] #160-160 MCP 3220V Mass 1 and 17
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-29', [r'scan_tof_2022-04-29-17-39-49'], False, 190, [0,1], True] #160-160 MCP 3800V Mass 2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-29', [r'scan_tof_2022-04-29-19-09-40'], False, 190, [0,1], True] #160-160 MCP 3800V Mass 2 and 16 long scan
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-30', [r'scan_tof_2022-04-30-16-49-32'], False, 190, [0,1], True] #120-120 MCP 3180V Mass 18
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-04-30', [r'scan_tof_2022-04-30-19-58-17'], False, 200, [0,1], True] #200-200 MCP 3800V Mass 2 and 16 long scan

    yearAndMonth = r'2022-10'
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-01', [r'scan_tof_2022-05-01-15-58-26'], False, 190, [0,1], True] #180-180 MCP 3800V Mass 2,16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-01', [r'scan_tof_2022-05-01-19-40-18'], False, 190, [0,1], True] #200-200 MCP 3800V Mass 2,16 !!!!!!!Good
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-01', [r'scan_tof_2022-05-01-22-48-45'], False, 190, [0,1], True] #180-180 MCP 3800V Mass 2,16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-02', [r'scan_tof_2022-05-02-11-16-58'], False, 190, [0,1], True] #100-100 MCP 3100V Mass18

    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-02', r'scan_tof_2022-05-02-13-46-03', False, 190, [0,1], True] #200-185 MCP 3800V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-02', [r'scan_tof_2022-05-02-16-19-42',r'scan_tof_2022-05-02-13-46-03'], False, 190, [0,1], True] #200-185 MCP 3100V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-02', [r'scan_tof_2022-05-02-20-04-38'], False, 190, [0,1], True] #200-200 MCP 3800V Mass2 and 16 !!!!!!!Good
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-03', [r'scan_tof_2022-05-03-10-15-09'], False, 190, [0,1], True] #160-160 MCP 3800V Mass2 and 16 !!!!!!!Good
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-03', [r'scan_tof_2022-05-03-12-47-32'], False, 190, [0,1], True] #160-160 MCP 3800V Mass2 and 16 !!!!!!!Good
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-03', [r'scan_tof_2022-05-03-15-35-46',r'scan_tof_2022-05-03-10-15-09'], False, 199, [2,4], True] #160-160 MCP 3800V Mass2 and 16 !!!!!!!Good
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-05', [r'scan_tof_2022-05-05-13-02-31'], False, 199, [0,1], True] #250-240 MCP 3800V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-05', [r'scan_tof_2022-05-05-14-27-54'], False, 199, [0,1], True] #220-200 MCP 3800V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-05', [r'scan_tof_2022-05-05-16-10-59'], False, 199, [0,1], True] #220-200 MCP 3800V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-05', [r'scan_tof_2022-05-05-17-25-16',r'scan_tof_2022-05-05-16-10-59'], False, 199, [0,1], True] #220-200 MCP 3800V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-05', [r'scan_tof_2022-05-05-19-03-54'], False, 199, [0,1], True] #180-165 MCP 3800V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-06', [r'scan_tof_2022-05-06-11-16-24'], False, 199, [0,1], True] #160-145 MCP 3800V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-06', [r'scan_tof_2022-05-06-14-59-14'], False, 199, [0,1], True] #160-145 MCP 3800V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-06', [r'scan_tof_2022-05-06-16-58-05 - Copy'], False, 199, [0,1], True] #140-220 MCP 3800V Mass2 and 16
    #date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-05-06', [r'scan_tof_2022-05-06-18-32-47'], False, 199, [0,1], True] #100-220 MCP 3800V Mass2 and 16

    date, filename, firstTry, sumNo, usefulRange, cu = [r'2022-10-10', [r'scan_tof_2022-10-10-21-27-42'], False, 20, [0,1], True] #100-220 MCP 3800V Mass2 and 16
    
    
    dataPath = pl.PurePath(dataPath)
    #if not os.path.exists(os.path.join(dataPath, yearAndMonth, date, filename)):
    #    os.mkdir(os.path.join(dataPath, yearAndMonth, date, filename))
    dataFile =  [os.path.join(dataPath, yearAndMonth, date, f+r'.hdf5') for f in filename]
    #massDictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_massDict'+r'.pkl')
    #interSpecDictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_interSpec'+r'.pkl')
    #interSpecDictFile2 = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_interSpec2'+r'.pkl')
    #FFTdictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_FFTdict'+r'.pkl')
    #SFTdictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_SFTdict'+r'.pkl')
    #delayFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_delay'+r'.pkl')
    #delayDictFile = os.path.join(dataPath, yearAndMonth, date, filename, filename+r'_delayDict'+r'.pkl')
    #files = [dataFile, massDictFile, interSpecDictFile, interSpecDictFile2, FFTdictFile, SFTdictFile, delayFile, delayDictFile]
    adelay = []
    for j in np.arange(1):
        #print(j)
        f = FFT_ionS(filename=dataFile, dcRange=1, massRange=[5,50])
        ff=[f]
    
        for x in ff:
            
            x.read_split(overNight = False, firstTry =firstTry, sumNo = sumNo, usefulRange = usefulRange, cu=cu)
            #x.read_multi()
            #check the spectrum
            #for i in range(0,4):
            #   plt.plot(x.spec[i],label=str(i))
            #   np.savetxt('spec.txt',x.spec[i])
            #plt.plot(x.spec[0],label=str(0))
            #plt.legend()
            #plt.show()
            plt.plot(np.arange(x.data.shape[0]),np.sum(x.data,1))
            #plt.plot(x.calculator.pixel2mass(np.arange(x.data.shape[0])),np.sum(x.data,1),label=None)
            #plt.xlabel('Atomic Mass (amu)')
            #plt.ylabel('log10(count)')
            #plt.plot(np.arange(sum.shape[0]),np.sum(sum[:,int(140/660*sum.shape[1]):],1),label='Mass[100,600]')
            #plt.legend()
            plt.show()
            
            x.mass_spectra()
 
 
            x.window(windowSize=150, direction='left',shift=j)
            #x.spectraBottle["H+"]=x.spectraBottle["H+"]/np.amax(x.spectraBottle["H+"][-500:])
            x.spectraBottle["O+"]=x.spectraBottle["O+"]/np.amax(x.spectraBottle["O+"][-500:])
            x.spectraBottle["H2+"]=x.spectraBottle["H2+"]/np.amax(x.spectraBottle["H2+"][-500:])
            #x.spectraBottle["O-H2"]=x.spectraBottle["O+"]-x.spectraBottle["H2+"]
            #x.spectraBottle["O+H2"]=x.spectraBottle["O+"]+x.spectraBottle["H2+"]
            #x.spectraBottle["H+H2"]=x.spectraBottle["H+"]+x.spectraBottle["H2+"]
           
            #x.useFilter(3400/33.35641*1e12,4000/33.35641*1e12)
 
            #x.rmvExp()
            #x.smooth(windowSize=11)
            x.show_Spectra()
            plt.show()
            #x.rebinS(factor=20)
            #x.STFTS(windowsize=int(200/600*3600))
            #x.smooth(windowSize=199)
            x.padding(paddingSize=33200)
            #x.useFilter(150/33*1e12,1600/33*1e12) 
            x.dataProcess()
            
            x.show_FFT()
            plt.show()
    #np.savetxt('phase_shift_160-160.dat',adelay)
    #plt.plot(adelay)
    #plt.show()
        
    