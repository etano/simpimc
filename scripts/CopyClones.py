import sys
import shutil

run1 = int(sys.argv[1])

try:
    move = int(sys.argv[2])
    files = sys.argv[3:]
except:
    move = 0
    files = sys.argv[2:]

for file0 in files:
    prefix = '.'.join(file0.split('.')[:-3])
    run0 = file0.split('.')[-3]
    clone = file0.split('.')[-2]
    file1 = prefix+'.'+str(run1)+'.'+clone+'.h5'
    print file1
    if move:
        shutil.move(file0,file1)
    else:
        shutil.copy(file0,file1)

