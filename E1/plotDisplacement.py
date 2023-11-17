import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('carbon_dioxide.csv',delimiter=',',skip_header=1)
freq = np.fft.fftfreq(len(array),0.001) #Should be twice the signal timestep
freq = np.fft.fftshift(freq)
signalSum = np.array(array[:, 1]+array[:, 2]+array[:, 3])


spectrum =  abs(np.fft.fft(array[:,2])) 
spectrum = np.fft.fftshift(spectrum)

pi = 3.141592

fig,ax = plt.subplots(3,1,figsize=(10,10))
ax[0].plot(array[:,0],array[:,1])
ax[0].plot(array[:,0],array[:,2])
ax[0].plot(array[:,0],array[:,3])
ax[0].set_ylabel('Displacement (Å)')
ax[0].set_xlabel('time (ps)')
ax[0].legend(["O","C", "O"])
ax[0].grid()

ax[1].plot(array[1:,0],array[1:,4])
ax[1].plot(array[1:,0],array[1:,5])
ax[1].plot(array[1:,0],array[1:,6])
ax[1].set_ylabel('Energy (eV)')
ax[1].set_xlabel('time (ps)')
ax[1].legend(["Ekin","Epot", "Etot"])
ax[1].grid()

ax[2].scatter(2565/33,0)
ax[2].scatter(1480/33,0)

ax[2].plot(freq,np.abs(spectrum))
ax[2].set_ylabel('Energy (eV)')
ax[2].set_xlabel('Frequency ')
ax[2].set_xlim(-150,150)
ax[2].set_ylim(0,30)
#ax[2].xaxis.set_ticks(np.arange(-100, 100, 10))

ax[2].grid()

fig.savefig('Co2.pdf')
