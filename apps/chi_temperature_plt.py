from numpy import load
from matplotlib import pyplot as plt

datas = load('chi_temperature.npy')
us = [1, 2.75, 3, 3.25, 4]
u_inds = [0,2,4]
colors = [plt.cm.jet(i/float(len(u_inds)-1)) for i in range(len(u_inds))]

fig = plt.figure()

ax = fig.add_subplot(1,2,1)
ax.set_title('local')
for u_ind, color in zip(u_inds, colors):
    data = datas[u_ind]
    u = us[u_ind]
    betas = data[0,:]
    temperatures = 1./betas
    chis_pm_loc = data[2,:]
    ax.plot(temperatures, chis_pm_loc, color = color)
ax.set_xlabel('$T$')
ax.set_xlim(0,.1)
#ax.set_ylim(.3,.45)
ax.set_ylabel('$\chi_{+-}(T)$')

ax = fig.add_subplot(1,2,2)
ax.set_title('nearest-neighbour')
for u_ind, color in zip(u_inds, colors):
    data = datas[u_ind]
    u = us[u_ind]
    betas = data[0,:]
    temperatures = 1./betas
    chis_pm_nn = data[3,:]
    ax.plot(temperatures, chis_pm_nn, label = str(u), color = color)
ax.set_xlabel('$T$')
ax.set_xlim(0,.1)
#ax.set_ylim(-.3,0)
plt.legend(loc = 'upper right', title = '$U$')

plt.tight_layout()
plt.savefig('chi_temperature.pdf', dpi=300)
plt.close()
