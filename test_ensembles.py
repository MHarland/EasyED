from ensembles import MicrocanonicalEnsemble, CanonicalEnsemble, GrandcanonicalEnsemble
from hamiltonians import Hubbard
from util import report

dimerHamilton = Hubbard([[0, -1], [-1, 0]], 1, verbose = True)
report('')
isolatedDimer = MicrocanonicalEnsemble(dimerHamilton)
isolatedDimer.calcOccupation()
report('n_micro: '+str(isolatedDimer.getTotalOccupation()))
report('')
heatedDimer = CanonicalEnsemble(dimerHamilton, 1)
heatedDimer.calcOccupation()
report('n_canonical: '+str(heatedDimer.getTotalOccupation()))
report('')
heatedStuffedDimer = GrandcanonicalEnsemble(dimerHamilton, 1, .5)
heatedStuffedDimer.calcOccupation()
report('n_grand: '+str(heatedStuffedDimer.getTotalOccupation()))
report('')
heatedStuffedDimer.setMuByFilling(2, 0, 2)
report('Numerical Mu for half-filling is '+str(heatedStuffedDimer.mu))
