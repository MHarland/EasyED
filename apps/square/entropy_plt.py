from numpy import load, log
from matplotlib import pyplot as plt

datas = load('entropy.npy')
us = [.75, 1]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot([0, 10], [log(2)]*2, color = 'gray', linestyle = 'dashed')
for data, u in zip(datas, us):
    betas = data[0,:]
    entropies = data[1,:]
    temperatures = 1./betas
    ax.plot(temperatures, entropies, label = str(u))
ax.set_xlabel('$T$')
ax.set_ylabel('$S(T)$')
ax.set_xlim(0,.1)
plt.legend()
plt.tight_layout()
plt.savefig('entropy.pdf', dpi=300)
plt.close()
