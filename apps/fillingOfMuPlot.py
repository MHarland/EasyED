from matplotlib import pyplot as plt, cm
from numpy import load

datas = load('fillingOfMu.npy')
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ulabels = ['0', '1', '2']
markers = ['o', 'x', '+']
for i, data, label, marker in zip(range(len(datas)), datas, ulabels, markers):
    ax.plot(data[0, :], data[1, :], color = cm.jet(i/float(len(datas)-1)), label = label)
ax.set_xlim([0, 4])
ax.set_ylim([0, 8])
ax.set_xlabel('$\mu$')
ax.set_ylabel('$N_{tot}$')
plt.legend()
plt.tight_layout()
fname = 'fillingOfMu.pdf'
plt.savefig(fname, dpi=300)
plt.close()
print fname, 'ready'
