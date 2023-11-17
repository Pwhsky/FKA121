import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('T10.csv',delimiter=',',skip_header=1)

fig,ax = plt.subplots(1,1,figsize=(15,10))
xlim = 1500
for i in range(1,32):
    plt.loglog(array[:,0],array[:,i],label = f"{i}")
ax.set_xlabel("time (ps)")
ax.set_ylabel("energy (Ev)")



fig.savefig('T10.pdf')
