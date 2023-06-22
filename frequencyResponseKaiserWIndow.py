import originpro as op
import os
my_path = os.path.abspath(__file__)
op.set_show(show=True)
windowFre = op.new_sheet('w',lname='windowFre')
from numpy.fft import fft, fftshift
import matplotlib.pyplot as plt
import numpy as np
window1 = np.kaiser(101, 10)
window2 = np.hanning(101)
window3 = np.hamming(101)
window4 = np.heaviside(np.arange(0,100,1), 1)
label=['Kaiser','Hann', 'Hamm', 'Step']
plt.plot(window1)

plt.title("Kaiser window")

plt.ylabel("Amplitude")

plt.xlabel("Sample")

plt.show()

plt.figure()

i=0
for window in [window1,window2,window3,window4]:
    A = fft(window, 2048) / 25.5

    mag = np.abs(fftshift(A))

    freq = np.linspace(-0.5, 0.5, len(A))

    response = 20 * np.log10(mag)

    response = np.clip(response, -100, 100)
    windowFre.from_list(0, freq, lname="Normalized frequency", axis='X')
    windowFre.from_list(i+1, response, lname=label[i], axis='Y')
    plt.plot(freq, response,label=label[i])
    
    i=i+1
plt.legend()
plt.title("Frequency response of Kaiser window")


plt.ylabel("Magnitude [dB]")


plt.xlabel("Normalized frequency [cycles per sample]")


plt.axis('tight')
(-0.5, 0.5, -100.0, ...) # may vary

plt.show()