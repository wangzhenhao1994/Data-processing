import h5py
import numpy as np
import pathlib as pl
import os
import datetime
from matplotlib import pyplot as plt
from matplotlib.pyplot import show
import originpro as op

op.set_show(show=True)


wks = op.new_sheet('w',lname='pulseSpectra')


now = datetime.datetime.now()
dataPath = pl.PureWindowsPath(
        #r'C:\Users\user\Desktop\Data_newTOF\2021-08\2021-08-16')#150-150
        r'C:\Users\user\Desktop\Data_newTOF\2022-06\2022-06-04')#170-170
        #r'C:\Users\user\Desktop\Data_newTOF\2021-06\2021-06-14')#70-70
        #r'C:\Users\user\Desktop\Data_newTOF\2021-06\2021-06-19')#70-50
        #r'C:\Users\user\Desktop\Data_newTOF\2021-06\2021-06-20')#70-40
datafile =  os.path.join(dataPath, r'spec.hdf5')#150-150
with h5py.File(datafile, 'r') as f:
    for key in f.keys():
        print(key)
