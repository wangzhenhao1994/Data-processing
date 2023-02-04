import numpy as np
import pathlib as pl
import os
from matplotlib import pyplot as plt
import qmsolve as qs
from scipy.interpolate import interp1d
import shutil


#interaction potential
path=pl.PureWindowsPath(r"D:\Program\simulation\javierelpianista-h2p-benchmarks-1cc95c4\discurves")
for file in os.listdir(path):
    r=np.loadtxt(os.path.abspath(os.path.join(path, file)))[:, 0]
    e=np.loadtxt(os.path.abspath(os.path.join(path, file)))[:, 1]
    if np.min(e)<e[0]:
        plt.plot(r,e)
        ff=os.path.splitext(file)[0]
        if not os.path.exists(os.path.join(path, ff)):
          os.makedirs(os.path.join(path,ff))
        np.savetxt(os.path.join(os.path.join(path, ff), 'pot_1.dat'),np.loadtxt(os.path.abspath(os.path.join(path, file)))[:, 0:2])
        shutil.copyfile(r'C:\Users\wangz\Downloads\qm_init.m', os.path.join(os.path.join(path, ff), r'qm_init.m'))
        shutil.copyfile(r'C:\Users\wangz\Downloads\qm_run.m', os.path.join(os.path.join(path, ff), r'qm_run.m'))

plt.show()





