import numpy as np
from matplotlib import pyplot as plt
import scipy
array = np.genfromtxt('radial.csv',delimiter=',',skip_header=1)


fig,ax = plt.subplots(figsize=(10,7))
#normalize = np.trapz(array[:,3],array[:,0])
min_g = min(array[40:-1,3])
r_index = np.where(array[:,3]==min_g)[0][0]
I_r = 4*np.pi* np.trapz(array[0:r_index,3]*(array[0:r_index,0]**2), array[0:r_index,0])
#print(r_index)
ax.scatter(array[r_index,0], min_g)
ax.plot(array[:,0],array[:,3],label="RDF")
ax.axvline( array[r_index,0])
#ax.plot(array[:,0],array[:,2],label="N_ideal")
#ax.plot(array[:,0],array[:,1],label="N_avg")

ax.set_xlabel('radius (Å) ',fontsize=15)
ax.set_ylabel('g(r)',fontsize=15)
ax.legend(loc="lower right")
fig.savefig('radial.pdf')
print("Correlation number 1 is: %f", I_r)
