import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('trajectory.dat')
plt.plot(data[:, 1], data[:, 2]) # Plot x versus y
plt.xlabel('x')
plt.xlabel('y')
plt.show()
