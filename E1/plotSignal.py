import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('signal.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots()
ax.plot(array[:,0],array[:,1])
ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()
fig.savefig('signal.pdf')