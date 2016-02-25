import sys
import os
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

# Get start cut
try:
    startCut = int(sys.argv[1])
    section = sys.argv[2]
    name = sys.argv[3]
    firstArg = 4
except:
    startCut = 0
    section = sys.argv[1]
    name = sys.argv[2]
    firstArg = 3
print section, name+'.eps'

# Plot traces
fig = plt.figure()
ax  = fig.add_subplot(111)

files = sys.argv[firstArg:]
count = np.array([0])
all_data = np.array([0])
for file in files:
    file_not_read = True
    failure_i = 0
    max_failure = 10
    while (file_not_read and failure_i < max_failure):
        try:
            f = h5.File(file,'r')
            data = np.array(f[section][startCut:])
            ax.plot(data)
            if data.size > all_data.size:
                all_data.resize(data.shape)
                count.resize(data.shape)
                all_data = all_data + data
                count = count + np.ones(data.shape)
            elif all_data.size > data.size:
                data2 = data.copy()
                data2.resize(all_data.shape)
                all_data = all_data + data2
                count2 = np.ones(data.shape)
                count2.resize(all_data.shape)
                count = count + count2
            else:
                all_data = all_data + data
                count = count + np.ones(data.shape)
            f.flush()
            f.close()
            file_not_read = False
        except IOError as e:
            print 'Trouble reading', section, 'in', file
            failure_i += 1
        except KeyError as e:
            print 'Data not found for', section, 'in', file
            file_not_read = False
    if (failure_i == max_failure):
        print 'Skipping', section, 'in', file
if count.size > 1 and all_data.size > 1:
    all_data /= count
    ax.plot(all_data,c='red',linewidth=5,label='Average')
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    # Show plot
    plt.title(section)
    plt.savefig(name+'.eps')
