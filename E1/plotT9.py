import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('T9.csv',delimiter=',',skip_header=1)

fig,ax = plt.subplots(1,1,figsize=(15,10))

for i in range(1,6):
    plt.loglog(array[:,0],array[:,i])

#plt.plot(array[:,0],array[:,1],label = "1")
#plt.plot(array[:,0],array[:,2],label = "2")
#plt.plot(array[:,0],array[:,3],label = "3")
#plt.plot(array[:,0],array[:,4],label = "4")
#plt.plot(array[:,0],array[:,5],label = "5")


plt.legend()



fig.savefig('T9.pdf')
