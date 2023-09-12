import numpy as np
#a=[1122,1285,1388,2244,3329]
#a=[523,813,1344,1743,2157,2778,3000,3464]
a=[526.94738,811.57589,884.65618,1146,1343.65171,1597.50957,1697.51418,2155.22759,2328.3125,2678.32864,3212.96868,3656.57888,1146.2067]
a = np.array(a,dtype='int')
dict = {}
for x in a:
    for y in a:
        dict[str(x)+'+'+str(y)] = x+y
        dict[str(x)+'-'+str(y)] = x-y
for key in dict.keys():
    for x in a:
        if np.abs(x-dict[key])<5:
            print(x,key)