from matplotlib import pyplot as plt
from numpy import array
from tetrahedron import tetrahedron

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

tetrahedron.calcG1([(('up',0),('up',0)), (('dn',0),('dn',0))])
tetrahedron.g1.setMesh(250, -1, 1)
sZ0 = .5 * (tetrahedron.g1.getRetarded((('up',0),('up',0))) - tetrahedron.g1.getRetarded((('dn',0),('dn',0))) )

ax.plot(tetrahedron.g1.getMesh(), sZ0.imag)
ax.set_xlabel('$\omega$')
ax.set_ylabel('$\mathrm{Im}\,S_{z,\,0}(\omega)$')
ax.set_ylim(min([-1, sZ0.imag.min()]), max([1, sZ0.imag.max()]))
plt.tight_layout()
plt.savefig('s_z.pdf', dpi=300)
plt.close()
