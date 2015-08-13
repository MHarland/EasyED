from itertools import product
from numpy import matrix, sum as nsum, zeros, argmin
from numpy.linalg import eigh
from operators import SingleParticleBasis, AnnihilationOperator
from util import diracDelta, dot

class Hamiltonian(SingleParticleBasis):
    def solve(self):
        self.eigenEnergies, self.eigenStates = eigh(self.h)

    def getGroundStateEnergy(self):
        return min(self.eigenEnergies)

    def getGroundStateAlgebraically(self):
        ind = argmin(self.eigenEnergies)
        return self.getStateAlgebraically(self.eigenStates[:, ind])

    def getGroundState(self):
        ind = argmin(self.eigenEnergies)
        return self.eigenStates[:, ind]

class Hubbard(Hamiltonian):
    def __init__(self, t, u):
        self.t = matrix(t)
        self.u = u
        self.spins = ['up', 'dn']
        self.sites = range(len(t))
        self.h = setHubbardMatrix(self.t, self.u, self.spins, self.sites)
        Hamiltonian.__init__(self, [self.spins, self.sites])

def setHubbardMatrix(t, u, spins, sites):
    c = AnnihilationOperator([spins, sites])
    ht = [t[i,j] * dot(c[s,i].H, c[s,j]) for s in spins for i,j in product(sites, sites)]
    hu = [.5 * u * diracDelta(s1, s2) * dot(c[s1,i].H, c[s1,i], c[s2,i].H, c[s2,i]) for i in sites for s1,s2 in product(spins, spins)]
    return nsum(ht + hu, axis = 0)
    
