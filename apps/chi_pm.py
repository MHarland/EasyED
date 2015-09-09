from matplotlib import pyplot as plt
from EasyED.operators import AnnihilationOperator
from EasyED.observables import DynamicObservable
from tetrahedron import tetrahedron

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax2 = fig.add_subplot(2,1,2)

c = AnnihilationOperator(tetrahedron.singleParticleBasis)

sPlus0 = c['up', 0].H.dot(c['dn', 0])
sMinus0 = c['dn', 0].H.dot(c['up', 0])
sMinus1 = c['dn', 1].H.dot(c['up', 1])
chi_pm = DynamicObservable({'loc': (sPlus0, sMinus0), 'nn': (sPlus0, sMinus1)}, True)
tetrahedron.setLehmannTermsDynamic(chi_pm)
chi_pm.setMesh(1000, -1, 1)

for ind in ['loc', 'nn']:
    ax1.plot(chi_pm.getMesh(), chi_pm.getCustom(ind, [1, -1], [1, 0]).imag, label = str(ind))
    #ax2.plot(chi_pm.getMesh(), chi_pm.getCustom(ind, [1, -1], [1, -1]).real, label = str(ind))

ax1.legend()
ax1.set_ylabel('$\mathrm{Im}\,\chi^{ret}_{+-}(\omega)$')
#ax2.set_ylabel('$\mathrm{Re}\,\chi_{+-}(\omega)$')
ax1.set_xlabel('$\omega$')
plt.tight_layout()
plt.savefig('chi_pm.pdf', dpi=300)
plt.close()
