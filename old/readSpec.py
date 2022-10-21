from matplotlib import pyplot as plt
import numpy as np
import h5py
import pathlib as pl

rootPath= pl.PureWindowsPath(r'C:\Users\user\Desktop\Data_newTOF')
filename=pl.PureWindowsPath(rootPath, r'2022-06\2022-06-06', r'scan_tof_2022-06-06-11-48-00.hdf5')
with h5py.File(filename, 'r+') as f:
    for key in f.keys():
        print(key)
    plt.plot(np.sum(f['data'],axis=0))
    plt.show()
