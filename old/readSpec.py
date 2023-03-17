from matplotlib import pyplot as plt
import numpy as np
import h5py
import pathlib as pl

rootPath= pl.PureWindowsPath(r'D:\DataProcessing')
filename=pl.PureWindowsPath(rootPath, r'2022-06\2022-06-01', r'spec.hdf5')
with h5py.File(filename, 'r+') as f:
    spec=0
    for key in f.keys():
        spec = np.array(f[key]) + spec
    plt.plot(spec)
    plt.show()
    np.savetxt('spec.txt',spec)
