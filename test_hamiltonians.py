from itertools import product
from hamiltonians import Hubbard
from matplotlib import pyplot as plt
from numpy import array, sqrt, sort
import numpy

from util import report

numpy.set_printoptions(suppress=True)
structure = Hubbard([[0, -1], [-1, 0]], 0)
#t = -1.
#structure = Hubbard([[0,t,t,t],[t,0,t,t],[t,t,0,t],[t,t,t,0]], 4, mu = 2)
#report(structure.matrix.toarray())
report(structure.blocksizes)
#print
#structure = Hubbard([[0, -1], [-1, 0]], 2, 1/sqrt(2) * array([[1,1],[1,-1]]))
#print structure.matrix.toarray()
#print
structure.solve()
#print structure.eigenEnergies
print sort(structure.eigenEnergies)
"""
for i in range(len(structure.eigenStates)):
    vstr = str()
    for j in range(len(structure.eigenStates)):
        val = structure.eigenStates[i,j]
        if abs(val) < .0001:
            vstr += ' 0.0'
        else:
            vstr += ' ' + str(val)[:4]
    print vstr
"""
print 'groundstateenergy: ', structure.getGroundStateEnergy()
print 'states:'
for state in structure.getGroundStatesAlgebraically():
    print state
    print
#print
#for state in structure.getGroundStates():
#
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
