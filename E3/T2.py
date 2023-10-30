import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('weighted.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots()
ax.scatter(array[:,0],array[:,1])
ax.set_xlabel('x ')
ax.set_ylabel('y')
ax.grid()
fig.savefig('T2.pdf')