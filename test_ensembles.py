from ensembles import MicrocanonicalEnsemble, CanonicalEnsemble, GrandcanonicalEnsemble
from hamiltonians import Hubbard

dimerHamilton = Hubbard([[0, -1], [-1, 0]], 1, verbose = True)
print
isolatedDimer = MicrocanonicalEnsemble(dimerHamilton)
isolatedDimer.calcOccupation()
print 'n_micro: ',isolatedDimer.getTotalOccupation()
print
heatedDimer = CanonicalEnsemble(dimerHamilton, 1)
heatedDimer.calcOccupation()
print 'n_canonical: ',heatedDimer.getTotalOccupation()
print
heatedStuffedDimer = GrandcanonicalEnsemble(dimerHamilton, 1, .5)
heatedStuffedDimer.calcOccupation()
print 'n_grand: ',heatedStuffedDimer.getTotalOccupation()
print
heatedStuffedDimer.setMu(2, 0, 2)
print 'Numerical Mu for half-filling is '+str(heatedStuffedDimer.mu)
