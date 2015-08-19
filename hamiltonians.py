from itertools import product
from numpy import array, sum as nsum, zeros, argmin
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from operators import SingleParticleBasis, AnnihilationOperator
from util import diracDelta, dot

class Hamiltonian(SingleParticleBasis):
    def __init__(self, singleParticleBasis, matrix):
        self.matrix = matrix
        SingleParticleBasis.__init__(self, singleParticleBasis)
        self.eigenEnergies = None
        self.eigenStates = None
        self.blocksizes = [self.fockspaceSize]
        self.sortN = None
        self.sortNSubspace()

    def sortNSubspace(self):
        self.blocksizes = list()
        self.sortN = list()
        nEigenvalues = [nsum([1 for digit in self.getOccupationRep(fockind) if digit == '1']) for fockind in range(self.fockspaceSize)]
        found = list()
        for fockspaceNr, n in enumerate(nEigenvalues):
            if n in found:
                self.sortN[found.index(n)].append(fockspaceNr)
            else:
                self.sortN.append(list())
                found.append(n)
                self.sortN[-1].append(fockspaceNr)
        self.blocksizes = [len(block) for block in self.sortN]
        self.sortN = [n for blockn in self.sortN for n in blockn]
        self.sortN = coo_matrix(([1]*len(self.sortN), (self.sortN, range(len(self.sortN)))), [self.fockspaceSize]*2)
        self.matrix = self.sortN.transpose().dot(self.matrix).dot(self.sortN)

    def sortNsSubspace(self):
        pass

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
    def __init__(self, t, u, siteSpaceTransformation = None, mu = 0):
        self.t = array(t) - identity(len(t)) * mu
        self.u = u
        self.spins = ['up', 'dn']
        self.sites = range(len(t))
        self.siteSpaceTransformation = siteSpaceTransformation
        hubbardMatrix = setHubbardMatrix(self.t, self.u, self.spins, self.sites, siteSpaceTransformation)
        Hamiltonian.__init__(self, [self.spins, self.sites], hubbardMatrix)

def setHubbardMatrix(t, u, spins, orbitals, siteSpaceTransformation):
    spins = range(len(spins))
    c = AnnihilationOperator([spins, orbitals])
    no = len(orbitals)
    ns = len(spins)
    spininds = range(len(spins))
    uMatrix = zeros([no,no,no,no,ns,ns])
    for i, j, k, l, s1, s2 in product(orbitals,orbitals,orbitals,orbitals,spins,spins):
        if i == k and j == l and i == j and s1 != s2:
            uMatrix[i, j, k, l, s1, s2] = u
    if siteSpaceTransformation != None:
        p = array(siteSpaceTransformation)
        t = p.transpose().dot(t).dot(p)
        temp = uMatrix.copy()
        for i, j, k, l, s1, s2 in product(orbitals,orbitals,orbitals,orbitals,spins,spins):
            uMatrix[i,j,k,l,s1,s2] = nsum([p[i,m] * p[j,n] * temp[m,n,o,q,s1,s2] * p.transpose()[o,l] * p.transpose()[q,k] for m,n,o,q in product(orbitals,orbitals,orbitals,orbitals)], axis = 0)
    ht = [t[i,j] * dot(c[s,i].transpose().dot(c[s,j])) for s in spins for i,j in product(orbitals, orbitals)]
    hu = [.5 * uMatrix[i, j, k, l, s1, s2] * c[s1,i].transpose().dot(c[s2, j].transpose()).dot(c[s2, l]).dot(c[s1, k]) for i,j,k,l in product(orbitals, orbitals, orbitals, orbitals) for s1, s2 in product(spins, spins) if s1 != s2]
    return nsum(ht + hu, axis = 0)
