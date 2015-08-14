from itertools import product
from numpy import matrix, sum as nsum, zeros, argmin
from numpy.linalg import eigh
from operators import SingleParticleBasis, AnnihilationOperator, SuperpositionState
from util import diracDelta, dot

class Hamiltonian(SingleParticleBasis):
    def __init__(self, singleParticleBasis, matrix):
        self.matrix = matrix
        SingleParticleBasis.__init__(self, singleParticleBasis)
        self.eigenEnergies = None
        self.eigenStates = None

    def solve(self):
        self.eigenEnergies, self.eigenStates = eigh(self.matrix)
    """
    def getGroundStateEnergy(self):
        return min(self.eigenEnergies)

    def getGroundStateAlgebraically(self):
        ind = argmin(self.eigenEnergies)
        psi0 = SuperpositionState(self.eigenStates[:, ind], self.singleParticleBasis)
        return psi0.getStateAlgebraically()

    def getGroundState(self):
        ind = argmin(self.eigenEnergies)
        psi0 = SuperpositionState(self.eigenStates[:, ind], self.singleParticleBasis)
        return psi0.getState()
    """
    def getSpectrum(self):
        degeneracies = list()
        energies = list()
        for e in self.eigenEnergies:
            if e in energies:
                degeneracies[energies.index(e)] += 1
            else:
                degeneracies.append(1)
                energies.append(e)
        return energies, degeneracies

class Hubbard(Hamiltonian):
    def __init__(self, t, u):
        self.t = matrix(t)
        self.u = u
        self.spins = ['up', 'dn']
        self.sites = range(len(t))
        hubbardMatrix = setHubbardMatrix(self.t, self.u, self.spins, self.sites)
        Hamiltonian.__init__(self, [self.spins, self.sites], hubbardMatrix)

def setHubbardMatrix(t, u, spins, sites):
    c = AnnihilationOperator([spins, sites])
    ht = [t[i,j] * dot(c[s,i].H, c[s,j]) for s in spins for i,j in product(sites, sites)]
    hu = [.5 * u * diracDelta(s1, s2) * dot(c[s1,i].H, c[s1,i], c[s2,i].H, c[s2,i]) for i in sites for s1,s2 in product(spins, spins)]
    return nsum(ht + hu, axis = 0)
    
