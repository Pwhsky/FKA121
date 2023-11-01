import numpy as np
import matplotlib.pyplot as plt
pi = 3.14159
N = 250
dt = 0.1
time = np.arange(0,250)*dt

signal = np.cos(time*2*pi*2 + pi/2)
#signal = np.cos(time*2*pi*2 + 0) + np.cos(2*pi*6*time)
fig,ax = plt.subplots()
ax.plot(time,signal)
ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()
fig.savefig('signal.pdf')
