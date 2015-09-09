from matplotlib import pyplot as plt
from numpy import load

datas = load('chi_temperature.npy')
betas = datas[0,:]
chis_pm_loc = datas[1,:]
chis_pm_nn = datas[2,:]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(1./betas, chis_pm_loc, label = 'loc')
ax.plot(1./betas, chis_pm_nn, label = 'nn')
ax.set_xlabel('$T$')
ax.set_ylabel('$\chi_{+-}(T)$')
ax.set_xlim(0,.3)
plt.legend()
plt.tight_layout()
plt.savefig('chi_temperature.pdf', dpi=300)
plt.close()
