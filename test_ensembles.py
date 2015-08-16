from ensembles import MicrocanonicalEnsemble, CanonicalEnsemble, GrandcanonicalEnsemble
from hamiltonians import Hubbard
from time import time, sleep

startTime = time()
dimerHamilton = Hubbard([[0, -1], [-1, 0]], 1)

isolatedDimer = MicrocanonicalEnsemble(dimerHamilton)
isolatedDimer.calcOccupation()
print 'n_micro: ',isolatedDimer.getTotalOccupation()

heatedDimer = CanonicalEnsemble(dimerHamilton, 1)
heatedDimer.calcOccupation()
print 'n_canonical: ',heatedDimer.getTotalOccupation()

heatedStuffedDimer = GrandcanonicalEnsemble(dimerHamilton, 1, .5)
heatedStuffedDimer.calcOccupation()
print 'n_grand: ',heatedStuffedDimer.getTotalOccupation()

duration = time() - startTime
sleep(1)
print 'took', duration, ' seconds'
