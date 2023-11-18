import numpy as np
from matplotlib import pyplot as plt
array = np.genfromtxt('task3.csv',delimiter=',',skip_header=1)


fig,ax = plt.subplots(figsize=(10,7))

ax.plot(array[:,0],array[:,2])
#
ax.set_xlabel('Timesteps ',fontsize=15)
ax.set_ylabel('Temperature (K)',fontsize=15)
ax.legend(loc="lower right")
fig.savefig('task3.pdf')