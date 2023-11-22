import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('block_averaging.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots()
ax.set_title("Statistical inefficency using block averaging")
ax.scatter(array[:,0],array[:,1])
ax.set_xlabel('block size ')
ax.set_ylabel('Statistical inefficency')
ax.grid()
fig.savefig('T7BlockAveraging.pdf')

array = np.genfromtxt('correlation.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots()
ax.set_title("Auto-Correlation function of MD simulation")
ax.plot(array[:,0],array[:,1])
ax.set_xlabel('time lag')
ax.set_ylabel('correlation')
ax.grid()
fig.savefig('T7Correlation.pdf')