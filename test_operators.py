from itertools import product
from numpy import sum as nsum, trace
from util import dot

from operators import SingleParticleBasis, AnnihilationOperator

inds = [['u', 'd'], [0, 1]]
basis =  SingleParticleBasis(inds)
print basis.getState(5)
print basis.getSingleParticleBasis()
print basis.getSingleParticleStateNr('d', 0)
print basis.getFockspaceNr(0,0,1,0)
print basis.getBinaryRep(4)


c = AnnihilationOperator([['up', 'dn'], range(2)])
print c.orderedSingleParticleStates
print c.blocksizes
print c.sortByN
print c['dn', 1]
#
#print c.element(('d', 'a'), 8, 12)
#print dot(c['u', 'a'].H, c['u', 'a'])
#n_tot = nsum([dot(c[s, i].H, c[s, i]) for s, i in product(['u', 'd'], ['a', 'b'])], axis = 0)
#print trace(dot(c['u', 'a'].H, c['u', 'a']))

