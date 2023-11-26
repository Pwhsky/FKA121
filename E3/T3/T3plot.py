import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('T3.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots(2,1,figsize=(15,7))

ax[0].plot(array[:,2],array[:,0])
ax[0].scatter(array[0,2],array[0,0],label="Starting point after burn")
ax[0].set_xlabel("time (microseconds)")
ax[0].set_ylabel("position (micrometers)")
ax[0].grid(True)
ax[0].legend()


ax[1].plot(array[:,2],array[:,1])
ax[1].scatter(array[0,2],array[0,1],label="Velocity after burn")
ax[1].set_xlabel("time (microseconds)")
ax[1].set_ylabel("velocity (micrometers/s)")
ax[1].grid(True)
ax[1].legend()
fig.suptitle("BD3 algorithm trajectory for high relaxation time, dt = 0.001")
fig.savefig("T3.pdf")