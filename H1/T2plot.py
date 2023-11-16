import numpy as np
from matplotlib import pyplot as plt
lowdt = np.genfromtxt('task2smalldt.csv',delimiter=',',skip_header=1)
highdt = np.genfromtxt('task2largedt.csv',delimiter=',',skip_header=1)

#bigdt = np.genfromtxt('task2bigdt.csv',delimiter=',',skip_header=1)

fig,ax = plt.subplots(figsize=(10,7))
temperature = lowdt[1:,1]*2/(3*8.6173*10**-5)
#ax.plot(array[:,0],array[:,1],label="kinetic energy")
#ax.plot(array[:,0],array[:,2],label="potential energy")
#mean1 = round(np.average(lowdt[:,4]))
#mean2 = round(np.average(highdt[:,4]))

ax.plot(range(1,100000),lowdt[1:,4],label= f"dt = 0.001")
#ax.plot(range(1,100000),highdt[1:,3],label= f"dt = 0.01")

#ax.plot(range(1,100),bigdt[1:,3],label= f"dt = 0.1")
#ax.plot(array[:,0],array[:,4],label="temperature")
ax.set_xlabel('Timesteps ',fontsize=15)
ax.set_ylabel('Temperature (K)',fontsize=15)
ax.legend(loc="lower right")
fig.savefig('task2temp.pdf')