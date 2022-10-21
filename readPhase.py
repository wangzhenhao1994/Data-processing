from scipy.signal import argrelextrema
import numpy as np

x = np.array([2, 1, 2, 3, 2, 0, 1, 0])
print(argrelextrema(x, np.greater))