import numpy as np
import matplotlib.pyplot as plt
pi = 3.14159
N = 255
dt = 0.1
time = np.arange(0,N)*dt
freq = np.fft.fftfreq(len(time),dt)
freq = np.fft.fftshift(freq)
#signal = np.cos(time*2*pi*2 + pi/2)
signal  = np.cos(time*2*pi*2 + 0) + np.cos(2*pi*6*time)
fourier_transform = np.fft.fft(signal)
fig,ax = plt.subplots()
ax.plot(freq,np.abs(fourier_transform))
ax.set_xlabel('freq (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()
fig.savefig('spectrum.pdf')
