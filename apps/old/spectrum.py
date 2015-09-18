from matplotlib import pyplot as plt
from tetrahedron import tetrahedron

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
es = list()
ds = list()
for e, d in zip(*tetrahedron.hamiltonian.getSpectrum()):
    es.append(e)
    ds.append(d)
    ax.plot([e]*2, [0, d], color = 'blue', zorder = 5)
xPlotRange = [min(es)-1,max(es)+1]
yPlotRange = [0,max(ds)]
for dPlotStep in range(max(ds)+1):
    ax.plot(xPlotRange, [dPlotStep]*2, color = 'gray', linestyle = '--', zorder = 1)
ax.set_xlim(*xPlotRange)
ax.set_ylim(*yPlotRange)
ax.set_xlabel('$\omega$')
ax.set_ylabel('$A(\omega)$')
plt.tight_layout()
plt.savefig('spectrum.pdf', dpi=300)
plt.close()