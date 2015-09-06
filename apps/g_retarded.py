from matplotlib import pyplot as plt
from numpy import array
from tetrahedron import tetrahedron

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
tetrahedron.g1.setMesh(250, 0, 1)

statePairs = [(('up', 0),('up', 0)), (('up', 0),('up', 1))]
for statePair in statePairs:
    tetrahedron.calcG1([statePair])
    tetrahedron.g1.setRetarded()
    ax.plot(tetrahedron.g1.getMesh(), tetrahedron.g1.getSpectralFunction(statePair), label = str(statePair))
ax.set_xlabel('$\omega$')
ax.set_ylabel('$-1/\pi\,\mathrm{Im}\,G^{ret}(\omega)$')
plt.legend()
plt.tight_layout()
plt.savefig('g_ret.pdf', dpi=300)
plt.close()
