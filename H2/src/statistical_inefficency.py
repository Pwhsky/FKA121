import numpy as np
import matplotlib.pyplot as plt



step = 3
P = np.genfromtxt('P.csv',delimiter=',')
U = np.genfromtxt('U.csv',delimiter=',')
T = np.genfromtxt('T.csv',delimiter=',')
r = np.genfromtxt('r.csv',delimiter=',')


Ct = np.genfromtxt('Ct.csv',delimiter=',')
Pt = np.genfromtxt('Pt.csv',delimiter=',')
Ut = np.genfromtxt('Ut.csv',delimiter=',')

sr = np.genfromtxt('sr.csv',delimiter=',')
sP = np.genfromtxt('sP.csv',delimiter=',')
sU = np.genfromtxt('sU.csv',delimiter=',')

data = [sU,sP,sr]
error = []
for i in range(0,3):

    #err = abs(estimated_autocorrelation(data[i]))
    err = abs(data[i])
    N = len(err)
    
    err = np.sqrt(err/N) *np.var(err)
    error.append(err)

fig,ax = plt.subplots(figsize=(12,7))

ax.plot(T,sU,label = "Internal energy",linestyle="--")

ax.plot(T,sP,label="Order parameter")





ax.plot(T,sr,label="Short range parameter",linestyle=":")



ax.set_xlabel('Temperature (K)')
ax.grid()
ax.legend()
fig.savefig('bars.pdf')
#ax[0, 0].errorbar(
#    T[::step],
#    U[::step],
#    yerr=error[0][::step],
#    label="Simulated",
#    fmt="o",
#    color="black",
#    capsize=3.5,  # Adjust the capsize for better visibility
#   elinewidth=2,  # Adjust the error bar line width  # Use a marker for data points
#    markersize=3,  # Increase marker size for better visibility
#    alpha=0.7,  # Add some transparency to reduce clutter
#)


#ax[1, 0].errorbar(
#    T[::step],
#    P[::step],
#    yerr=error[1][::step],
#    label="Simulated",
#   fmt="o",
#    color="black",
#    capsize=3.5,  # Adjust the capsize for better visibility
#    elinewidth=2,  # Adjust the error bar line width  # Use a marker for data points
#    markersize=3,  # Increase marker size for better visibility  
#    alpha=0.7,  
#    # Add some transparency to reduce clutter
#)


#ax[1, 1].errorbar(
#    T[::step],
#    r[::step],
#   yerr=error[2][::step],
#    label="Simulated",
#   fmt="o",
#    color="black",
#    capsize=3.5,  # Adjust the capsize for better visibility
#    elinewidth=2,  # Adjust the error bar line width  # Use a marker for data points
#   markersize=3,  # Increase marker size for better visibility
#   alpha=0.7,  # Add some transparency to reduce clutter
#)
#fig.savefig('bars.pdf')