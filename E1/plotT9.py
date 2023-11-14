import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('T9.csv',delimiter=',',skip_header=1)

fig,ax = plt.subplots(1,1,figsize=(15,10))
xlim = 1500
plt.plot(array[:xlim,0],array[:xlim,1],label = "1")
plt.plot(array[:xlim,0],array[:xlim,2],label = "2")
plt.plot(array[:xlim,0],array[:xlim,3],label = "3")
plt.plot(array[:xlim,0],array[:xlim,4],label = "4")
plt.plot(array[:xlim,0],array[:xlim,5],label = "5")
plt.legend()



fig.savefig('T9.pdf')
