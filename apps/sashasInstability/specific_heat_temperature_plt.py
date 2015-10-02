from numpy import load
from matplotlib import pyplot as plt
#from matplotlib_to_latex import set_poster_parameters as set_mpl

datas = load('specific_heat_temperature.npy')
mus = [.26,.27,.28]
inds = [0,1,2]
colors = [plt.cm.jet(i/float(max(len(inds), 2)-1)) for i in range(len(inds))]

#set_mpl()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for ind, color in zip(inds, colors):
    mu = mus[ind]
    data = datas[ind]
    betas = data[0,:]
    temperatures = 1./betas
    specificHeat = data[1,:]
    ax.plot(temperatures, specificHeat, color = color, label = str(mu))
ax.set_xlabel('$T$')
ax.set_xlim(0,.1)
ax.set_ylabel('$C(T)$')
plt.legend(loc = 'upper right', title = '$\mu$')
plt.tight_layout()
plt.savefig('specific_heat_temperature.pdf', dpi=300)
plt.close()
