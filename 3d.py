import numpy as np
import matplotlib.pyplot as plt
import originpro as op
op.set_show(show=True)
wks = op.new_sheet('w',lname='time')
i=0
j=0
phi = [1.45,1.71,1,0.83,3.16,3.44]
amp = (np.array([2.5608,2.68038,2.68751,2.48369,2.44122,3.4])-2.4)*0.5
for f in [521,805,1334,2145,3201,3644]:
    t = np.linspace(0,100,1300)
    y = np.zeros(1300)+f
    z = np.sin(2*np.pi/33356.40952*f*t+phi[j])*amp[j]+amp[j]/2
    j=j+1
    wks.from_list(i, t, axis='X')
    wks.from_list(i+1, y, axis='Y')
    wks.from_list(i+2, z, axis='Z')
    i=i+3