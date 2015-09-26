from itertools import product
from hamiltonians import Hubbard
from matplotlib import pyplot as plt
from numpy import array, sqrt, sort
import numpy

from util import report

numpy.set_printoptions(suppress=True)
structure = Hubbard([[0, -1], [-1, 0]], 1)

report(structure.blocksizes)
structure.solve()
print structure.eigenEnergies

print 'groundstateenergy: ', structure.getGroundStateEnergy()
print 'groundstates:'
for state in structure.getGroundStatesAlgebraically():
    print state
    print

print 'states:'
for energyGroup, energy in zip(structure.getStatesEnergySortedAlgebraically(), structure.getEnergies()):
    print 'E = ', energy
    print
    for state in energyGroup:
        print state
        print
    print
    print

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
es = list()
for e, d in zip(*structure.getSpectrum()):
    es.append(e)
    ax.plot([e]*2, [0, d], color = 'blue')
ax.set_xlim(min(es)-1,max(es)+1)
plt.tight_layout()
plt.savefig('testSpectrum.pdf', dpi=300)
plt.close()
