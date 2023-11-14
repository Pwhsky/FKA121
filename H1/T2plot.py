import numpy as np
from matplotlib import pyplot as plt
array = np.genfromtxt('task2.csv',delimiter=',',skip_header=1)


fig,ax = plt.subplots()
#ax.plot(array[:,0],array[:,1],label="kinetic energy")
#ax.plot(array[:,0],array[:,2],label="potential energy")
ax.plot(array[:,0],array[:,3],label= "total energy")
#ax.plot(array[:,0],array[:,4],label="temperature")
ax.set_xlabel('time ($ps$)',fontsize=14)
ax.set_ylabel('Energy (eV)',fontsize=14)
ax.legend()
fig.savefig('task2.pdf')