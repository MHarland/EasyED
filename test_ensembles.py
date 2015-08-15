from ensembles import CanonicalEnsemble, GrandcanonicalEnsemble
from hamiltonians import Hubbard
from time import time

startTime = time()
dimerHamilton = Hubbard([[0, -1], [-1, 0]], 1)
#heatedDimer = CanonicalEnsemble(dimerHamilton, 10)
#print heatedDimer.occupation('up', 0)
heatedStuffedDimer = GrandcanonicalEnsemble(dimerHamilton, 10, -5)
print heatedStuffedDimer.totalOccupation()
heatedStuffedDimer = GrandcanonicalEnsemble(dimerHamilton, 10, 15)
print heatedStuffedDimer.totalOccupation()
duration = time() - startTime
print 'took', duration
