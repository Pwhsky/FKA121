import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('weighted.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots()
ax.scatter(array[:,1],array[:,0])
ax.set_xlabel('x ')
ax.set_ylabel('y')
ax.grid()
fig.savefig('T2.pdf')