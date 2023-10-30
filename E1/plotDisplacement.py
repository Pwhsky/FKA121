import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('coupled_oscillator.csv',delimiter=',',skip_header=1)
spectrum = np.genfromtxt('powerspectrum.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots(3,1,figsize=(10,10))
ax[0].plot(array[:,0],array[:,1])
ax[0].plot(array[:,0],array[:,2])
ax[0].plot(array[:,0],array[:,3])
ax[0].set_ylabel('Displacement (Å)')
ax[0].set_xlabel('time (ps)')
ax[0].legend(["O","C", "O"])
ax[0].grid()

ax[1].plot(array[:,0],array[:,4])
ax[1].plot(array[:,0],array[:,5])
ax[1].plot(array[:,0],array[:,6])
ax[1].set_ylabel('Energy (eV)')
ax[1].set_xlabel('time (ps)')
ax[1].legend(["Ekin","Epot", "Etot"])
ax[1].grid()

spectrumSum = np.array(spectrum[:, 0]+spectrum[:, 1]+spectrum[:, 2])
spectrum = np.array(spectrum)*33.3564095
ax[2].plot(spectrum[:,3],spectrumSum )
ax[2].set_ylabel('Energy (eV)')
ax[2].set_xlabel('wavenumber (cm^-1)')
ax[2].set_xlim(-5000,5000)

ax[2].grid()



fig.savefig('coupled_oscillator.pdf')
