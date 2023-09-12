import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
#amp  =1.967e+05
#n    =       1 
#phi  =   0
#w0   =    3600
#dw=4.244736138859549
#x=np.linspace(3550,3650,300)
#
#for deltaPhi in np.linspace(0,2*np.pi,40):
#    model = amp*np.sinc(np.pi*(x-w0)/dw)*np.exp(-1j*(np.pi*(x-w0)/dw+deltaPhi))
#    plt.plot(x,np.imag(model),'b')
#    plt.plot(x,np.real(model),'r')
#    plt.plot(x,np.abs(model),'g')
#    plt.show()
x=np.linspace(0,1300,13000)
for deltaPhi in np.linspace(0,2*np.pi,40):
    y = np.cos(x/9.12125-deltaPhi)
    y2=y*np.hanning(np.size(y))
    y2=np.concatenate((y2,np.zeros(130000)))
    y=np.concatenate((y,np.zeros(130000)))
    fft = np.fft.rfft(y)
    plt.plot(np.real(fft),'r')
    plt.plot(np.imag(fft))
    plt.title(str(deltaPhi))
    plt.plot(np.abs(fft),'g')
    plt.xlim([200,300])
    plt.show()
#
    #fft2 = np.fft.rfft(y2)
    #ref=np.concatenate((np.hanning(np.size(y)),np.zeros(130000)))
    #refP = np.fft.rfft(ref)
    ##plt.plot(np.real(fft2),'r')
    ##plt.plot(np.imag(fft2))
    ##plt.title(str(deltaPhi))
    #plt.plot(np.angle(fft2)-np.angle(refP))
    #plt.xlim([100,300])
    #plt.plot(np.abs(fft2),'g')
    #plt.show()