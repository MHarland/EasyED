from numpy import loadtxt

from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load, save, array, where, linspace

almostZero = 10**(-10)
mu = .27
beta = 20
results = list()
u = 3
t = -1
r = .3
h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
tetrahedron.calcEigensystem()
energies = list(tetrahedron.hamiltonian.eigenEnergies)

energiesSasha = loadtxt('/afs/physnet.uni-hamburg.de/users/th1_li/group-th1_li/alichten/VERTEX/ED_Xi_Om_4s/TEST/U3N3/energy.dat')[:,1]

for j, eS in enumerate(energiesSasha):
    found = False
    for i, e in enumerate(energies):
        if abs(eS - e) < almostZero:
            found = True
            del energies[i]
            break
    report('found nr '+str(j))
    assert found, str(eS)+' not found'
