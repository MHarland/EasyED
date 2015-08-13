from hamiltonians import Hubbard

dimer = Hubbard([[0, -1], [-1, 0]], 1.5)
dimer.solve()
print dimer.getGroundStateEnergy()
print
print dimer.getGroundStateAlgebraically()
