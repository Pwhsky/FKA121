import numpy as np
import matplotlib.pyplot as plt

def estimated_autocorrelation(x):
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result

c_avg = np.genfromtxt(f'C.csv', delimiter=',') 
a = np.max(c_avg)
print(a)
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
    
    err = np.sqrt(err/N) * np.sqrt(np.var(err))
    error.append(err)

fig,ax = plt.subplots(2,2,figsize=(12,7))

ax[0,0].plot(T,U,label = "Simulated")

ax[0,0].plot(Ut[:,0],Ut[:,1],label = "Mean field theory")
ax[0,0].set_ylabel('Internal energy (eV)')



ax[0,1].plot(T[1:],c_avg,label="Simulated")

ax[0,1].plot(Ct[:,0],Ct[:,1],label = "Mean field theory")
ax[0,1].set_ylabel('Heat capacity (eV/K)')


ax[1,0].plot(T,P,label="Simulated")
ax[1,0].plot(Pt[:,0],Pt[:,1],label = "Mean field theory")
ax[1,0].set_ylabel('Order Parameter')




ax[1,1].plot(T,r,label="Simulated")
ax[1,1].set_ylabel('short range parameter')


for i in range(2):
    for j in range(2):

        ax[i,j].axvline(x=745, color='red', linestyle='--', label='T_c = 745 K')
        ax[i,j].legend()
        ax[i,j].set_xlabel('Temperature (K)')
        ax[i,j].grid()
        ax[i,j].legend()

fig.savefig('T2.pdf')
#ax[0, 0].errorbar(
#    T[::step],
#    U[::step],
#    yerr=error[0][::step],
#    label="Simulated",
#    fmt="o",
#    color="black",
#    capsize=3.5,  # Adjust the capsize for better visibility
#    elinewidth=2,  # Adjust the error bar line width  # Use a marker for data points
#    markersize=3,  # Increase marker size for better visibility
#    alpha=0.7,  # Add some transparency to reduce clutter
#)


#ax[1, 0].errorbar(
#    T[::step],
#    P[::step],
#    yerr=error[2][::step]/75,
#    label="Simulated",
#    fmt="o",
#    color="black",
#    capsize=3.5,  # Adjust the capsize for better visibility
#    elinewidth=2,  # Adjust the error bar line width  # Use a marker for data points
#    markersize=3,  # Increase marker size for better visibility
#    alpha=0.7,  # Add some transparency to reduce clutter
#)


#ax[1, 1].errorbar(
#    T[::step],
#    r[::step],
 #   yerr=error[3][::step]/75,
  #  label="Simulated",
 #   fmt="o",
 #   color="black",
 #   capsize=3.5,  # Adjust the capsize for better visibility
  #  elinewidth=2,  # Adjust the error bar line width  # Use a marker for data points
 #   markersize=3,  # Increase marker size for better visibility
 #   alpha=0.7,  # Add some transparency to reduce clutter
#)
