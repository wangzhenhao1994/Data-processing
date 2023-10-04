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
    y=np.concatenate((y,np.zeros(130000)))
    y2=-y
    fft = np.fft.rfft(y)
    fft2 = np.fft.rfft(y2)

    plt.plot(np.real(fft),'r')
    plt.plot(np.imag(fft),'r--')


    
    plt.plot(np.real(fft2),'g')
    plt.plot(np.imag(fft2),'g--')
    plt.title(str(deltaPhi))

    plt.xlim([100,300])
    #plt.plot(np.abs(fft2),'g')
    plt.show()
#x=np.linspace(0,10,13000)
#y = np.cos((x/4)*2*np.pi-0.5*np.pi)
#y2 = np.cos((x/4)*2*np.pi-0.5*np.pi)
#plt.figure(figsize=(8, 2), dpi=600)
#ax = plt.subplot(111)
#ax.plot(x, y, 'b', linewidth = 2)
#ax.plot(x, y2*0.5,'r',linewidth = 2)
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#ax.spines['left'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
#plt.yticks([])
#plt.xticks([])
#plt.savefig('test.png',transparent = True, dpi=600)
#plt.show()