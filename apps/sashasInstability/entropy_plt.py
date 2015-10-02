from numpy import load, log
from matplotlib import pyplot as plt
from matplotlib_to_latex import set_poster_parameters as set_mpl

datas = load('entropy.npy')
mus = [.26,.27,.28]
inds = [0,1,2]
colors = [plt.cm.jet(i/float(len(inds)-1)) for i in range(len(inds))]

set_mpl()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#ax.plot([0, 10], [log(2)]*2, color = 'gray', linestyle = 'dashed')
#ax.plot([0, 10], [4*log(2)]*2, color = 'gray', linestyle = 'dotted')
for ind, color in zip(inds, colors):
    data = datas[ind]
    mu = mus[ind]
    betas = data[0,:]
    entropies = data[1,:]
    temperatures = 1./betas
    ax.plot(temperatures, entropies/log(4), label = str(mu), color = color)
ax.set_xlabel('$T$')
ax.set_ylabel('$S(T)/\mathrm{ln}(4)$')
#ax.set_xlim(0,1)
#ax.set_ylim(0,8)
plt.legend(title = '$\mu$', loc = 'lower right')
plt.tight_layout()
plt.savefig('entropy.pdf', dpi=300)
plt.close()
