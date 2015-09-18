from numpy import load
from matplotlib import pyplot as plt
from matplotlib_to_latex import set_poster_parameters as set_mpl

datas = load('specific_heat_temperature.npy')
us = [1, 2.25, 2.5, 2.75, 3, 3.25, 4]
u_inds = [3,4,5]
u_inds = [3,4,5]
colors = [plt.cm.jet(i/float(max(len(u_inds), 2)-1)) for i in range(len(u_inds))]

set_mpl()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for u_ind, color in zip(u_inds, colors):
    u = us[u_ind]
    data = datas[u_ind]
    betas = data[0,:]
    temperatures = 1./betas
    specificHeat = data[1,:]
    ax.plot(temperatures, specificHeat, color = color, label = str(u))
ax.set_xlabel('$T$')
ax.set_xlim(0,.1)
ax.set_ylim(0,400)
ax.set_ylabel('$C(T)$')
plt.legend(loc = 'upper right', title = '$U$')
plt.tight_layout()
plt.savefig('specific_heat_temperature.pdf', dpi=300)
plt.close()
