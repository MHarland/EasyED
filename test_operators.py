from itertools import product
from operators import AnnihilationOperator
from numpy import sum as nsum, trace
from util import dot

c = AnnihilationOperator([['u', 'd'], ['a', 'b']])
print c.orderedSingleParticleStates
#print c['d',1,'p']
#print c.spsNr('u', 3, 's')
print c.fsNr(1,1,1,1)
print c.element(('d', 'a'), 8, 12)
print dot(c['u', 'a'].H, c['u', 'a'])
n_tot = nsum([dot(c[s, i].H, c[s, i]) for s, i in product(['u', 'd'], ['a', 'b'])], axis = 0)
print trace(dot(c['u', 'a'].H, c['u', 'a']))
