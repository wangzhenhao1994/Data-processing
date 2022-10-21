import os
directory = os.path.abspath(r'C:\Users\user\Desktop\Data_newTOF')
for root, dirs, files in os.walk(directory):
    #print(root)
    path = root.split(os.sep)
    for file in files:
        if '.txt' in file:
            with open(os.path.join(root, file), 'r') as f:
                a = f.read()
                if r'Reflected beam: 200mW' in a:
                    if r'Transmitted beam: 200mW' in a:
                        print(file)