import numpy as np
from matplotlib import pyplot as plt
array = np.genfromtxt('task4.csv',delimiter=',',skip_header=1)


fig,ax = plt.subplots(figsize=(10,7))

ax.plot(array[:,0],array[:,4],label="Atom 1 position")
ax.plot(array[:,0],array[:,5],label="Atom 2 position")

ax.set_xlabel('Time (ps) ',fontsize=15)
ax.set_ylabel('Position (Å)',fontsize=15)
ax.legend(loc="lower right")
fig.savefig('task4Atoms.pdf')


fig,ax = plt.subplots(figsize=(10,7))
ax.plot(array[:,0],array[:,1])

ax.set_xlabel('Time (ps) ',fontsize=15)
ax.set_ylabel('Temperature (K)',fontsize=15)
fig.savefig('task4temperature.pdf')

fig,ax = plt.subplots(figsize=(10,7))
ax.plot(array[:,0],array[:,2])

ax.set_xlabel('Time (ps) ',fontsize=15)
ax.set_ylabel('Pressure (bar)',fontsize=15)
fig.savefig('task4pressure.pdf')

fig,ax = plt.subplots(figsize=(10,7))
ax.plot(array[:,0],array[:,3])

ax.set_xlabel('Time (ps) ',fontsize=15)
ax.set_ylabel('Lattice constant (Å)',fontsize=15)
fig.savefig('task4lattice.pdf')