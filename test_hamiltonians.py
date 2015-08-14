from hamiltonians import Hubbard
from matplotlib import pyplot as plt
from numpy import trace, sum as nsum, matrix
from operators import AnnihilationOperator
from util import dot

dimer = Hubbard([[0, -1], [-1, 0]], 1)
dimer.solve()
#print dimer.getGroundStateEnergy()
#print dimer.getGroundStateAlgebraically()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.scatter(*dimer.getSpectrum())
plt.savefig('spectrum.pdf')
plt.close()
