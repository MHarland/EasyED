from matplotlib import pyplot as plt
from EasyED.operators import AnnihilationOperator
from EasyED.observables import DynamicObservable, lehmannSumDynamic
from tetrahedron import tetrahedron


fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax2 = fig.add_subplot(2,1,2)

c = AnnihilationOperator(tetrahedron.singleParticleBasis)
chi_pm = DynamicObservable(True)
chi_pm.partitionFunction = tetrahedron.getPartitionFunction()
sPlus0 = c['up', 0].H.dot(c['dn', 0])
sMinus0 = c['dn', 0].H.dot(c['up', 0])
sMinus1 = c['dn', 1].H.dot(c['up', 1])
noms, denoms = tetrahedron.getLehmannTermsDynamic(sPlus0, sMinus0)
chi_pm.lehmannNominators.update({'loc': noms})
chi_pm.lehmannDenominators.update({'loc': denoms})
noms, denoms = tetrahedron.getLehmannTermsDynamic(sPlus0, sMinus1)
chi_pm.lehmannNominators.update({'nn': noms})
chi_pm.lehmannDenominators.update({'nn': denoms})
chi_pm.setMesh(250, 0, 1)

for ind in ['loc', 'nn']:
    ax1.plot(chi_pm.getMesh(), chi_pm.getCustom(ind, [1, -1], [1, -1]).imag, label = str(ind))
    #ax2.plot(chi_pm.getMesh(), chi_pm.getCustom(ind, [1, -1], [1, -1]).real, label = str(ind))

ax1.legend()
ax1.set_ylabel('$\mathrm{Im}\,\chi^{c}_{+-}(\omega)$')
#ax2.set_ylabel('$\mathrm{Re}\,\chi_{+-}(\omega)$')
ax1.set_xlabel('$\omega$')
plt.tight_layout()
plt.savefig('chi_pm.pdf', dpi=300)
plt.close()
