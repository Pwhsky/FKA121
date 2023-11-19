import numpy as np
from matplotlib import pyplot as plt
array = np.genfromtxt('task3.csv',delimiter=',',skip_header=1)


fig,ax = plt.subplots(figsize=(10,7))

ax.plot(array[:,0],array[:,4],label="atom 1")
ax.plot(array[:,0],array[:,5],label="atom 2")
ax.plot(array[:,0],array[:,6],label="atom 3")

ax.set_xlabel('Time (ps) ',fontsize=15)
ax.set_ylabel('Position (Å)',fontsize=15)
ax.legend(loc="lower right")
fig.savefig('task3Atoms.pdf')


fig,ax = plt.subplots(figsize=(10,7))
ax.plot(array[:,0],array[:,1])

ax.set_xlabel('Time (ps) ',fontsize=15)
ax.set_ylabel('Temperature (K)',fontsize=15)
fig.savefig('task3temperature.pdf')

fig,ax = plt.subplots(figsize=(10,7))
ax.plot(array[:,0],array[:,2])

ax.set_xlabel('Time (ps) ',fontsize=15)
ax.set_ylabel('Pressure(eV/$Å^3$)',fontsize=15)
fig.savefig('task3pressure.pdf')

fig,ax = plt.subplots(figsize=(10,7))
ax.plot(array[:,0],array[:,3])

ax.set_xlabel('Time (ps) ',fontsize=15)
ax.set_ylabel('Lattice constant (Å)',fontsize=15)
fig.savefig('task4lattice.pdf')