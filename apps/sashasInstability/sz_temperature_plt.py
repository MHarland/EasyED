from numpy import load
from matplotlib import pyplot as plt

datas = load('sz_temperature.npy')
us = [3]
u_inds = [0]
colors = [plt.cm.jet(i/float(max(len(u_inds)-1, 1))) for i in range(len(u_inds))]

fig = plt.figure()

ax = fig.add_subplot(1,1,1)
for u_ind, color in zip(u_inds, colors):
    data = datas[u_ind]
    u = us[u_ind]
    betas = data[0,:]
    temperatures = 1./betas
    sz_loc = data[2,:]
    ax.plot(temperatures, sz_loc, color = color)
ax.set_xlabel('$T$')
ax.set_xlim(0,.1)
ax.set_ylim(-.1,.1)
#ax.set_ylim(.3,.45)
ax.set_ylabel('$<S_{z}>$')
plt.tight_layout()
plt.savefig('sz_temperature.pdf', dpi=300)
plt.close()
