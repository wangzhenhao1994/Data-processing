import numpy as np 
import matplotlib.pyplot as plt 

n = 5000
T = 2*np.pi
dw = 2 * np.pi / T

t = np.linspace(0, 1000*np.pi, n, endpoint=False)
x = np.sin(t)
fftx = np.fft.rfft(x)
freq = np.fft.rfftfreq(n) * n * dw

amps = np.abs(fftx) * 2 /  n
angs = np.angle(fftx)
#angs[amps < 1] = 0


_, ax = plt.subplots(3, 1, sharex=True)
ax[0].plot(t, x)
ax[1].plot(freq, amps)
ax[2].plot(freq, angs)
plt.show()