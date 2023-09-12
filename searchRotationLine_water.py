import numpy

for i in range(0,20,2):
    for j in range(20):
        for k in range(20):
            if abs(27.8770027*i+14.51200*j+9.28500*k-(812-527))<3:
                print(i,j,k)
 	 	