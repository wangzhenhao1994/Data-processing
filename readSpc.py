import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from pyspectra.readers.read_spc import read_spc
spc=read_spc(r'D:\Users\wangz\Desktop\water_spc\WTERK95.SPC')
spc.plot()
plt.xlabel("nm")
plt.ylabel("Abs")
plt.grid(True)
print(spc.head())
plt.show()
