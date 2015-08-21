from hamiltonians import Hubbard
from matplotlib import pyplot as plt
from numpy import array, sqrt, sort
import numpy

from util import report

numpy.set_printoptions(suppress=True)
structure = Hubbard([[0, -1], [-1, 0]], 4, mu = 2)
#t = -1.
#structure = Hubbard([[0,t,t,t],[t,0,t,t],[t,t,0,t],[t,t,t,0]], 4, mu = 2)
#report(structure.matrix.toarray())
#report(structure.blocksizes)
#print
#structure = Hubbard([[0, -1], [-1, 0]], 2, 1/sqrt(2) * array([[1,1],[1,-1]]))
#print structure.matrix.toarray()
#print
structure.solve()
print sort(structure.eigenEnergies)
#print structure.eigenStates
print 'groundstateenergy: ', structure.getGroundStateEnergy()
print 'states:'
for state in structure.getGroundStatesAlgebraically():
    print state
#print
#for state in structure.getGroundStates():
#


















    print state


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for e, d in zip(*structure.getSpectrum()):
    ax.plot([e]*2, [0, d], color = 'blue')
plt.savefig('testSpectrum.pdf', dpi=300)
plt.close()
