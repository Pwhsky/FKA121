import numpy as np
import matplotlib.pyplot as plt

C = np.genfromtxt('C.csv',delimiter=',',skip_header=1)
P = np.genfromtxt('P.csv',delimiter=',',skip_header=1)
U = np.genfromtxt('U.csv',delimiter=',',skip_header=1)

fig,ax = plt.subplots(3,1,figsize=(7,10))
ax[0].plot(U[:,0],U[:,1],label = "Internal energy")
ax[0].set_ylabel('Internal energy (eV)')
ax[0].set_xlabel('Temperature (K)')


step = 1
ax[1].plot(C[::step,0],C[::step,1],label="Heat capacity")

ax[1].set_ylabel('Heat capacity (eV/K)')
ax[1].set_xlabel('Temperature (K)')


ax[2].plot(P[:,0],P[:,1],label="Order parameter")
ax[2].set_ylabel('Order Parameter')
ax[2].set_xlabel('Temperature (K)')


for x in ax:
    x.axvline(x=905.15, color='red', linestyle='--', label='T_c = 905.15 K')
    x.legend()

fig.savefig('T1.pdf')
