from operators import AnnihilationOperator
from numpy import dot

c = AnnihilationOperator([['u', 'd'], ['a', 'b']])
print c.orderedSingleParticleStates
#print c['d',1,'p']
#print c.spsNr('u', 3, 's')
print c.fsNr(1,1,1,1)
print c.element(('d', 'a'), 8, 12)
print dot(c['u', 'a'].H, c['u', 'a'])
