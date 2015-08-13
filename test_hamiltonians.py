from hamiltonians import Hubbard
from matplotlib import pyplot as plt

dimer = Hubbard([[0, -1], [-1, 0]], 1)
dimer.solve()
print dimer.getGroundStateEnergy()
print
print dimer.getGroundStateAlgebraically()
print
print dimer.getGroundState()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.scatter(*dimer.getSpectrum())
plt.savefig('spectrum.pdf')
plt.close()
