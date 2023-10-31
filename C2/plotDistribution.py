import numpy as np
import matplotlib.pyplot as plt

array = np.genfromtxt('distribution.csv',delimiter=',',skip_header=1)
fig,ax = plt.subplots()
ax.scatter(array[:],array[:])
fig.savefig('distribution.pdf')