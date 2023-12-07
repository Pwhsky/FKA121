import numpy as np
from matplotlib import pyplot as plt
array = np.genfromtxt('radial.csv',delimiter=',',skip_header=1)


fig,ax = plt.subplots(figsize=(10,7))
normalize = np.trapz(array[:,3],array[:,0])

ax.axhline(y = 115, color = 'r', linestyle = '--') 
ax.plot(array[:,0],array[:,3],label="RDF")
#ax.plot(array[:,0],array[:,2],label="N_ideal")
#ax.plot(array[:,0],array[:,1],label="N_avg")

print(np.trapz(array[:,3],array[:,0]/normalize))
ax.set_xlabel('Distance (Å) ',fontsize=15)
ax.set_ylabel('g(r) ',fontsize=15)
ax.legend(loc="lower right")
fig.savefig('radial.pdf')
