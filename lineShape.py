import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

def fitFunc(x,phi):
    return np.sinc(np.pi*x)*np.exp(1j*(np.pi*x+phi))

x=np.linspace(-2,2,1000)
plt.subplot(1,2,1)
plt.plot(x,-np.real(fitFunc(x,-0.7*np.pi)),label='Real')
plt.plot(x,-np.imag(fitFunc(x,-0.7*np.pi)),label='Imag')
plt.legend()
plt.subplot(1,2,2)
plt.plot(x,np.real(fitFunc(x,-0.7*np.pi)),label='Real')
plt.plot(x,np.imag(fitFunc(x,-0.7*np.pi)),label='Imag')
plt.legend()
plt.show()