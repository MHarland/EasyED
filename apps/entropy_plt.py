from numpy import load, log
from matplotlib import pyplot as plt
from matplotlib_to_latex import set_poster_parameters as set_mpl

datas = load('entropyb.npy')
us = [1,2.25,2.5,2.75,3,3.25,4]
us = [1,2,3,4]
u_inds = [0,1,2,3]
colors = [plt.cm.jet(i/float(len(u_inds)-1)) for i in range(len(u_inds))]

set_mpl()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#ax.plot([0, 10], [log(2)]*2, color = 'gray', linestyle = 'dashed')
#ax.plot([0, 10], [4*log(2)]*2, color = 'gray', linestyle = 'dotted')
for u_ind, color in zip(u_inds, colors):
    data = datas[u_ind]
    u = us[u_ind]
    betas = data[0,:]
    entropies = data[2,:]
    temperatures = 1./betas
    ax.plot(temperatures, entropies/log(2), label = str(u), color = color)
ax.set_xlabel('$T$')
ax.set_ylabel('$S(T)/\mathrm{ln}(2)$')
ax.set_xlim(0,3)
ax.set_ylim(0,8)
plt.legend(title = '$U$', loc = 'lower right')
plt.tight_layout()
plt.savefig('entropyb.pdf', dpi=300)
plt.close()
