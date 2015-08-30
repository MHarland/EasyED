from matplotlib import pyplot as plt
from ensembles import GrandcanonicalEnsemble
from hamiltonians import Hubbard

dimerHamilton = Hubbard([[0, -1], [-1, 0]], 4, verbose = True)
print
dimer = GrandcanonicalEnsemble(dimerHamilton, 1, 2)
statePair = (('up', 1),('up', 1))
dimer.calcG1([statePair])

dimer.g1.setMesh(500, -6, 6)
dimer.g1.setRetarded()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(dimer.g1.getMesh(), dimer.g1.getSpectralFunction(statePair))
plt.tight_layout()
plt.savefig('testRetarded.pdf', dpi=300)
plt.close()

dimer.g1.setMatsubaraMesh(20, 1)
dimer.g1.setMatsubara()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(dimer.g1.getMatsubaraMesh(), dimer.g1.getMatsubara(statePair).real, label = 'Re')
ax.plot(dimer.g1.getMatsubaraMesh(), dimer.g1.getMatsubara(statePair).imag, label = 'Im')
plt.legend()
plt.tight_layout()
plt.savefig('testMatsubara.pdf', dpi=300)
plt.close()
