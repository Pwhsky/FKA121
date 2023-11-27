import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('T3.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots(figsize=(15,7))
ax.plot(array[:,2],array[:,0])
ax.scatter(array[0,2],array[0,0],label="Starting point after burn")
ax.set_xlabel("time (microseconds)")
ax.set_ylabel("position (micrometers)")
ax.grid(True)
ax.legend()
fig.savefig('T3position.pdf')

fig,ax = plt.subplots(figsize=(15,7))
ax.plot(array[:,2],array[:,1])
ax.scatter(array[0,2],array[0,1],label="Velocity after burn")
ax.set_xlabel("time (microseconds)")
ax.set_ylabel("velocity (micrometers/s)")
ax.grid(True)
ax.legend()
fig.savefig('T3velocity.pdf')