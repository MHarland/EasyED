from matplotlib import pyplot as plt
from EasyED.operators import AnnihilationOperator
from EasyED.observables import DynamicObservable
from tetrahedron import tetrahedron


fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)

c = AnnihilationOperator(tetrahedron.singleParticleBasis)
chi_zz = DynamicObservable(True)
chi_zz.partitionFunction = tetrahedron.getPartitionFunction()
sZ0 = .5 * (c['up', 0].H.dot(c['up', 0]) - c['dn', 0].H.dot(c['dn', 0]))
sZ1 = .5 * (c['up', 1].H.dot(c['up', 1]) - c['dn', 1].H.dot(c['dn', 1]))
noms, denoms = tetrahedron.getLehmannTermsDynamic(sZ0, sZ0)
chi_zz.lehmannNominators.update({'loc': noms})
chi_zz.lehmannDenominators.update({'loc': denoms})
noms, denoms = tetrahedron.getLehmannTermsDynamic(sZ0, sZ1)
chi_zz.lehmannNominators.update({'nn': noms})
chi_zz.lehmannDenominators.update({'nn': denoms})
chi_zz.setMesh(1000, 0, .25)

for ind in ['loc', 'nn']:
    ax1.plot(chi_zz.getMesh(), chi_zz.getCustom(ind, [1, 1], [1, 1]).imag, label = str(ind))
    ax2.plot(chi_zz.getMesh(), chi_zz.getCustom(ind, [1, 1], [1, 1]).real, label = str(ind))


ax1.legend()
ax1.set_ylabel('$\mathrm{Im}\,\chi_{zz}(\omega)$')
ax2.set_ylabel('$\mathrm{Re}\,\chi_{zz}(\omega)$')
ax2.set_xlabel('$\omega$')
plt.tight_layout()
plt.savefig('chi_zz.pdf', dpi=300)
plt.close()
