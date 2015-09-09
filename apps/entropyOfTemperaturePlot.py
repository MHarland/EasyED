from matplotlib import pyplot as plt
from numpy import load, log, exp

datas = load('entropyOfTemperature.npy')
betas = datas[0,:]
entropies = datas[1,:]
mus = datas[4,:]
x = 1./betas

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot([0, .3], [log(2)]*2, color = 'gray', linestyle = 'dashed')
ax.plot(x, entropies)
ax.set_xlabel('$T$')
ax.set_ylabel('$S(T)$')
#ax.set_ylim(0,4.5)
ax.set_xlim(0,.3)
plt.tight_layout()
plt.savefig('entropyOfTemperature.pdf', dpi=300)
plt.close()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(x, mus)
ax.set_xlabel('$T$')
ax.set_ylabel('$\mu(T)$')
ax.set_xlim(0,.3)
plt.tight_layout()
plt.savefig('muOfTemperature.pdf', dpi=300)
plt.close()
