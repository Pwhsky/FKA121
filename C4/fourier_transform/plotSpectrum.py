import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('powerspectrum.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots()
ax.plot(array[:,1],array[:,0])
ax.set_xlabel('freq (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()
fig.savefig('powerspectrum.pdf')
