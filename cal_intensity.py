import numpy as np
from decimal import Decimal
def cal_intensity(power, diameter, duration):
    repetition = 3000         
    power = power*1e-3 #power in W                       
    duration = duration*1e-15
    #diameter = 8                             #diameter of the laser spot in pixel
    diameter = diameter * 4.6 * 1e-4         #diameter in centimeter, 4.6 um is the pixel size
    area = np.pi*diameter*diameter/4
    intensity = power/repetition/area/duration
    return intensity
if __name__ == '__main__':
    print('%.1E' % Decimal(cal_intensity(250,9.2,5)))
    print('%.1E' % Decimal(cal_intensity(200,10.5,5)))
    print('%.1E' % Decimal(cal_intensity(230,9.7,5)))