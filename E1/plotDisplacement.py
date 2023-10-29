import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('coupled_oscillator.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots()
ax.plot(array[:,0],array[:,1])
ax.plot(array[:,0],array[:,2])
ax.plot(array[:,0],array[:,3])
ax.set_ylabel('Displacement (Å)')
ax.set_xlabel('time (ps)')
ax.grid()
fig.savefig('coupled_oscillator.pdf')
