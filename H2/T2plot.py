import numpy as np
import matplotlib.pyplot as plt


c_avg = np.zeros(0)  # Initialize an empty array

for i in range(10):
    new_c = np.genfromtxt(f'C{i}.csv', delimiter=',')
    c_avg = np.append(c_avg, new_c)

# Reshape c_avg to a 2D array
c_avg = c_avg.reshape(10, -1)

# Compute the mean element-wise along axis 0
mean_array = np.mean(c_avg, axis=0)

P = np.genfromtxt('P.csv',delimiter=',')
U = np.genfromtxt('U.csv',delimiter=',')
T = np.genfromtxt('T.csv',delimiter=',')
r = np.genfromtxt('r.csv',delimiter=',')

Ct = np.genfromtxt('Ct.csv',delimiter=',')
Pt = np.genfromtxt('Pt.csv',delimiter=',')
Ut = np.genfromtxt('Ut.csv',delimiter=',')




fig,ax = plt.subplots(2,2,figsize=(10,7))

ax[0,0].plot(T,U,label = "Simulated")
ax[0,0].plot(Ut[:,0],Ut[:,1],label = "Theory")
ax[0,0].set_ylabel('Internal energy (eV)')
ax[0,0].set_xlabel('Temperature (K)')
ax[0,0].grid()
ax[0,0].legend()

ax[0,1].plot(T[1:],mean_array,label="Simulated")
ax[0,1].plot(Ct[:,0],Ct[:,1],label = "Theory")
ax[0,1].set_ylabel('Heat capacity (eV/K)')
ax[0,1].set_xlabel('Temperature (K)')
ax[0,1].grid()
ax[0,1].legend()

ax[1,0].plot(T,P,label="Simulated")
ax[1,0].plot(Pt[:,0],Pt[:,1],label = "Theory")
ax[1,0].set_ylabel('Order Parameter')
ax[1,0].set_xlabel('Temperature (K)')
ax[1,0].grid()
ax[1,0].legend()


ax[1,1].plot(T,r,label="Simulated")
ax[1,1].set_ylabel('short range parameter')
ax[1,1].set_xlabel('Temperature (K)')
ax[1,1].grid()
ax[1,1].legend()

fig.savefig('T2.pdf')
