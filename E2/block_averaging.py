#%%
import numpy as np
import matplotlib.pyplot as plt


#%% 
data = np.genfromtxt('MC.txt', delimiter=',', skip_header=1)

M = len(data)                       # number of data points
B_array = np.arange(50, 3000, 50)   # size of each block
n_s = np.zeros(len(B_array))

var_f = np.sum(data**2)/M - (np.sum(data)/M)**2

for n , B in enumerate(B_array):
    print("Block size: " + str(B))
    Mb = int(M/B)               # number of blocks
    F_blocks = np.zeros(Mb)

    for j in range(Mb):
        for i in range(B):
            F_blocks[j] += data[j*B+i]
        F_blocks[j] = F_blocks[j]/B

    # calculate  variance of F_blocks (var[F])
    var_F = np.sum(F_blocks**2)/Mb - (np.sum(F_blocks)/Mb)**2
    n_s[n] = B*var_F/var_f


fig, ax = plt.subplots()
ax.plot(B_array, n_s)
plt.show()
