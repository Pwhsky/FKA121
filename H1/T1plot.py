import numpy as np
from matplotlib import pyplot as plt
array = np.genfromtxt('task1.csv',delimiter=',',skip_header=1)
array[:,1] = array[:,1]**3

model = np.poly1d(np.polyfit(array[:,1], array[:,0], 3))
polyline = np.linspace(min(array[:,1]), max(array[:,1]), 50)
fitted_values = model(polyline)
fig,ax = plt.subplots()
ax.scatter(array[:,1],array[:,0],label="Computed point")
ax.plot(polyline,fitted_values,label="Quadratic fit",linestyle="dashed")

ax.set_xlabel('Unit cell volume ($Å^3$)',fontsize=14)
ax.set_ylabel('Energy (eV/unit cell)',fontsize=14)
ax.legend()
fig.savefig('task1.pdf')