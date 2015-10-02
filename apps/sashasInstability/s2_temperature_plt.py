from numpy import load
from matplotlib import pyplot as plt

datas = load('s2_temperature.npy')
us = [3]
u_inds = [0]
colors = [plt.cm.jet(i/float(max(len(u_inds)-1, 1))) for i in range(len(u_inds))]

fig = plt.figure()

ax = fig.add_subplot(1,2,1)
ax.set_title('local')
for u_ind, color in zip(u_inds, colors):
    data = datas[u_ind]
    u = us[u_ind]
    betas = data[0,:]
    temperatures = 1./betas
    s2_loc = data[2,:]
    ax.plot(temperatures, s2_loc, color = color)
ax.set_xlabel('$T$')
ax.set_xlim(0,.1)
ax.set_ylim(.4,.6)
ax.set_ylabel('$<S^{2}>$')

ax = fig.add_subplot(1,2,2)
ax.set_title('nearest-neighbour')
for u_ind, color in zip(u_inds, colors):
    data = datas[u_ind]
    u = us[u_ind]
    betas = data[0,:]
    temperatures = 1./betas
    s2_nn = data[3,:]
    ax.plot(temperatures, s2_nn, label = str(u), color = color)
ax.set_xlabel('$T$')
ax.set_xlim(0,.1)
ax.set_ylim(-.2,0)
plt.legend(loc = 'upper right', title = '$U$')

plt.tight_layout()
plt.savefig('s2_temperature.pdf', dpi=300)
plt.close()
