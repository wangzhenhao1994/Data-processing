#import matlab.engine
#eng = matlab.engine.start_matlab()
#import matlab
import numpy as np
from matplotlib import pyplot as plt
from read_data_osc_daq import FFT_ionS
from decimal import Decimal
import pathlib as pl
from math import ceil
import os
from PyWrapOrigin.PyWrapOrigin import PyWrapOrigin
class timeFreAnalysis(FFT_ionS):
    def __init__(self) -> None:
        self.molecule = 'H2O'
        self.mass = [1,17,18]
        self.intensity = [1.9E+14, 2.3E+14, 8.0E+14]
        self.stage = 'piezo'
        self.rootPath= pl.PureWindowsPath(r'C:\Users\user\Desktop\Data_newTOF\plotdata')
        self.timeTrace = {}

        self.scanTime = 5
        self.trueScanTime = ceil(self.scanTime/1.6)*1.6
        self.scanLength = self.scanTime*1000
        self.delay = np.arange(self.scanLength)/self.scanLength*((self.scanTime)/self.trueScanTime)*99*2*2*3.33564*10**-15

    def readTimeTrace(self):
        for mass in self.mass:
            for i in self.intensity:
                label = 'Mass '+str(mass)
                self.timeTrace[self.molecule+str(mass)+'_'+str(i)]=np.loadtxt(pl.PureWindowsPath(self.rootPath, str('Time')+label+self.molecule+self.stage+str('%.1E' % Decimal(i))+str('.dat')))

    def cpsd(self, x, y, fs):
        __len=min(np.size(x),np.size(y))
        x=x.flatten()
        y=y.flatten()
        x=x[-__len:]
        y=y[-__len:]
        print(np.shape(x))
        x=matlab.double(x.tolist(),is_complex=False)
        y=matlab.double(y.tolist(),is_complex=False)
        fs=matlab.double([fs])
        emp=matlab.double([])
        Cxy=np.array(eng.mscohere(x,y,emp,emp,emp,fs)).flatten()
        plt.plot(Cxy)
        plt.show()
        Pxy=np.array(eng.cpsd(x,y,emp,emp,emp,fs)).flatten()
        plt.plot(Pxy)
        plt.show()

if __name__ == '__main__':
    tfa=timeFreAnalysis()
    tfa.readTimeTrace()
    plt.plot(tfa.timeTrace[tfa.molecule+str(1)+'_'+str(tfa.intensity[0])])
    plt.show()
    #tfa.cpsd(tfa.timeTrace[tfa.molecule+str(1)+'_'+str(tfa.intensity[0])],tfa.timeTrace[tfa.molecule+str(18)+'_'+str(tfa.intensity[0])],fs=1/(tfa.delay[1]-tfa.delay[0]))
    #for mass in tfa.mass:
    #        for i in tfa.intensity:
    #            label = 'Mass '+str(mass)
    #            plt.plot(tfa.timeTrace[tfa.molecule+str(mass)+'_'+str(i)])
    #            plt.show()


