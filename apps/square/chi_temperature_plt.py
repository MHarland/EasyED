from numpy import load
from matplotlib import pyplot as plt

datas = load('chi_temperature.npy')
us = [0.75, 1, 2, 4]
colors = [plt.cm.jet(i/float(len(us)-1)) for i in range(len(us))]

fig = plt.figure()

ax = fig.add_subplot(1,2,1)
ax.set_title('local')
for data, u, color in zip(datas, us, colors):
    betas = data[0,:]
    temperatures = 1./betas
    chis_pm_loc = data[1,:]
    ax.plot(temperatures, chis_pm_loc, color = color)
ax.set_xlabel('$T$')
ax.set_xlim(0,.1)
ax.set_ylim(.3,.45)
ax.set_ylabel('$\chi_{+-}(T)$')

ax = fig.add_subplot(1,2,2)
ax.set_title('neartest-neighbour')
for data, u, color in zip(datas, us, colors):
    betas = data[0,:]
    temperatures = 1./betas
    chis_pm_nn = data[2,:]
    ax.plot(temperatures, chis_pm_nn, label = str(u), color = color)
ax.set_xlabel('$T$')
ax.set_xlim(0,.1)
ax.set_ylim(-.3,0)
plt.legend(loc = 'lower right')

plt.tight_layout()
plt.savefig('chi_temperature.pdf', dpi=300)
plt.close()
