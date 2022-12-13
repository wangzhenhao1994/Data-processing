import os
import re
directory = os.path.abspath(r'D:\DataProcessing')
for root, dirs, files in os.walk(directory):
    #print(root)
    path = root.split(os.sep)
    for file in files:
        if '.txt' in file:
            with open(os.path.join(root, file), 'r') as f:
                a = f.read()
                if r'Ch10 is the reference of Ch2, Ch4 and Ch6' in a:
                    print(file.replace('.txt',''))
                    if r'overlap' in a:
                        print(int(re.findall(r'\d+', a[0:45])[1]),int(re.findall(r'\d+', a[0:45])[0]))
                        try:
                            print(float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)", a[45:70])[1]),float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)", a[45:70])[0]))
                        except IndexError:
                            print(float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)", a[45:70])[0]))
                        print('\n')
                    else:
                        print(int(re.findall(r'\d+', a[0:45])[0]),int(re.findall(r'\d+', a[0:45])[1]))
                        try:
                            print(float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)", a[45:70])[0]),float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)", a[45:70])[1]))
                        except IndexError:
                            print(float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)", a[45:70])[0]))
                        print('\n')