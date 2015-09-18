from numpy import load, log
from matplotlib import pyplot as plt
from matplotlib_to_latex import set_poster_parameters as set_mpl

datas = load('entropy.npy')
us = [1,2.25,2.5,2.75,3,3.25,4]
u_inds = [0,4,6]
colors = [plt.cm.jet(i/float(len(u_inds)-1)) for i in range(len(u_inds))]

set_mpl()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot([0, 10], [log(2)]*2, color = 'gray', linestyle = 'dashed')
for u_ind, color in zip(u_inds, colors):
    data = datas[u_ind]
    u = us[u_ind]
    betas = data[0,:]
    entropies = data[2,:]
    temperatures = 1./betas
    ax.plot(temperatures, entropies, label = str(u), color = color)
ax.set_xlabel('$T$')
ax.set_ylabel('$S(T)$')
ax.set_xlim(0,.1)
#ax.set_ylim(0.5,.9)
plt.legend(title = '$U$', loc = 'upper left')
plt.tight_layout()
plt.savefig('entropy.pdf', dpi=300)
plt.close()
