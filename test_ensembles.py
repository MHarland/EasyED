from ensembles import CanonicalEnsemble, GrandcanonicalEnsemble
from hamiltonians import Hubbard

dimerHamilton = Hubbard([[0, -1], [-1, 0]], 1)
heatedDimer = CanonicalEnsemble(dimerHamilton, 10)
print heatedDimer.occupation('up', 0)
heatedStuffedDimer = GrandcanonicalEnsemble(dimerHamilton, 100, -5)
print heatedDimer.occupation('up', 0)
