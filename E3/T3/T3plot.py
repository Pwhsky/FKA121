import numpy as np
import matplotlib.pyplot as plt

ltht = np.genfromtxt('01.csv',delimiter=',',skip_header=1)
ltlt = np.genfromtxt('00.csv',delimiter=',',skip_header=1)
htlt = np.genfromtxt('10.csv',delimiter=',',skip_header=1)
htht = np.genfromtxt('11.csv',delimiter=',',skip_header=1)
time = htht[:,2]

data = [ltlt,htlt,ltht,ltlt]


fig,ax = plt.subplots(2,2,figsize=(10,7))
ax[0,0].plot(time,ltlt[:,0],label="99.8 kPa, dt = 0.001")
ax[0,1].plot(time,htlt[:,0],label="2.75 kPa, dt = 0.001")
ax[1,0].plot(time,ltht[:,0],label="99.8 kPa, dt = 0.005")
ax[1,1].plot(time,htht[:,0],label="2.75 kPa, dt = 0.005")

for i in range(2):
    for j in range(2):
        ax[i,j].set_xlabel("Time (ms)")
        ax[i,j].set_ylabel("position (micrometers)")
        ax[i,j].grid()
        ax[i,j].legend()


fig.savefig('T3position.pdf')

fig,ax = plt.subplots(2,2,figsize=(10,7))
ax[0,0].plot(time,ltlt[:,1],label="99.8 kPa, dt = 0.001")
ax[0,1].plot(time,htlt[:,1],label="2.75 kPa, dt = 0.001")
ax[1,0].plot(time,ltht[:,1],label="99.8 kPa, dt = 0.005")
ax[1,1].plot(time,htht[:,1],label="2.75 kPa, dt = 0.005")

for i in range(2):
    for j in range(2):
        ax[i,j].set_xlabel("Time (ms)")
        ax[i,j].set_ylabel("Velocity (micrometers/s)")
        ax[i,j].grid()
        ax[i,j].legend()


fig.savefig('T3velocity.pdf')

