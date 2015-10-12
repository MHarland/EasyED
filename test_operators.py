from itertools import product
from numpy import sum as nsum, trace
from util import dot

from operators import SingleParticleBasis, AnnihilationOperator

inds = [['u', 'd'], [0, 1]]
basis =  SingleParticleBasis(inds)
print basis.getSingleParticleBasis()
print basis.orderedSingleParticleStates
print 'u0 d0 d1:'
print basis.getStateAlgebraically(13)
print basis.getFockspaceNr((1,0,1,1))
print basis.getOccupationRep(13)
print 'up1:'
print basis.getFockspaceNr((0,1,0,0))
print basis.getSingleParticleStateNr('u', 1)
print basis.getOccupationRep(2)
print 'all single particle states:'
for sps in basis.orderedSingleParticleStates:
    print basis.getFockspaceNr(singleParticleState = sps)
print

c = AnnihilationOperator([['u', 'd'], range(2)])
print c.orderedSingleParticleStates
print c['u', 1].toarray()
print
print c['u', 1].transpose().dot(c['u', 1]).toarray()
print nsum([c[s, i].transpose().dot(c[s, i]) for s, i in product(['u', 'd'], range(2))], axis = 0)

