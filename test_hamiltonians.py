from hamiltonians import Hubbard
from matplotlib import pyplot as plt
from numpy import array, sqrt
import numpy

numpy.set_printoptions(suppress=True)
dimer = Hubbard([[0, -1], [-1, 0]], 1)
print dimer.matrix.toarray()
print 'blocks:', dimer.blocksizes
print
dimer = Hubbard([[0, -1], [-1, 0]], 0, 1/sqrt(2) * array([[1,1],[1,-1]]))
print dimer.matrix.toarray()
print
dimer = Hubbard([[0, -1], [-1, 0]], 2)
print dimer.matrix.toarray()
print
dimer.solve()
#print dimer.getGroundStateEnergy()
#print dimer.getGroundStateAlgebraically()

#fig = plt.figure()
#ax = fig.add_subplot(1from hamiltonians import Hubbard
#plt.close()
