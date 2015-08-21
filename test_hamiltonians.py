from hamiltonians import Hubbard
from matplotlib import pyplot as plt
from numpy import array, sqrt
import numpy

from util import report

numpy.set_printoptions(suppress=True)
#structure = Hubbard([[0, -1], [-1, 0]], 2)
t = -1.
structure = Hubbard([[0,t,t,t],[t,0,t,t],[t,t,0,t],[t,t,t,0]], 2)
report(structure.matrix.toarray())
#report(structure.blocksizes)
#print
#structure = Hubbard([[0, -1], [-1, 0]], 2, 1/sqrt(2) * array([[1,1],[1,-1]]))
#print structure.matrix.toarray()
#print
structure.solve()
#report(structure.eigenEnergies)
#report(structure.eigenStates)
#print structure.getGroundStateEnergy()
#print structure.getGroundStateAlgebraically()

#fig = plt.figure()
#ax = fig.add_subplot(1from hamiltonians import Hubbard
#plt.close()
