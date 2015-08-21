from itertools import product, izip
from numpy import array, sum as nsum, zeros, argmin, identity, asmatrix
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from time import time

from blocks import BlockMatrix
from operators import SingleParticleBasis, AnnihilationOperator
from util import scatter_list, allgather_list, report

class Hamiltonian(SingleParticleBasis):
    def __init__(self, singleParticleBasis, matrix, verbose = False):
        self.matrix = matrix
        SingleParticleBasis.__init__(self, singleParticleBasis)
        self.eigenEnergies = None
        self.eigenStates = None
        self.blocksizes = [self.fockspaceSize]
        self.sortN = None
        self.fockBasis = range(self.fockspaceSize)
        self.sortNSubspace()
        self.verbose = verbose

    def sortNSubspace(self):
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
        self.fockBasis = [n for blockn in self.sortN for n in blockn]
        self.sortN = coo_matrix(([1]*self.fockspaceSize, (range(self.fockspaceSize), self.fockBasis)), [self.fockspaceSize]*2)
        self.matrix = self.sortN.dot(self.matrix).dot(self.sortN.transpose())

    def solve(self):
        report('Solving the Hamiltonian...', self.verbose)
        t0 = time()
        self.eigenEnergies = list()
        self.eigenStates = list()
        hBlocks = BlockMatrix(self.blocksizes)
        rows, cols = self.matrix.nonzero()
        for val, i, j in izip(self.matrix.data, rows, cols):
            hBlocks[i, j] = val
        hBlocks_scat = scatter_list(hBlocks.datablocks)
        e_scat = list()
        v_scat = list()
        for block in hBlocks_scat:
            e, v = eigh(block)
            e_scat += list(e)
            v_scat += list(v)
        self.eigenEnergies = allgather_list(e_scat)
        v = allgather_list(v_scat)
        self.eigenStates = embedV(v, self.blocksizes)
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)
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
    def __init__(self, t, u, siteSpaceTransformation = None, mu = 0, verbose = False):
        self.verbose = verbose
        self.sortNs = list()
        self.t = array(t) - identity(len(t)) * mu
        self.u = u
        self.spins = ['up', 'dn']
        self.sites = range(len(t))
        self.siteSpaceTransformation = siteSpaceTransformation
        report('Setting up the Hubbard Hamiltonian...', self.verbose)
        t0 = time()
        hubbardMatrix = setHubbardMatrix(self.t, self.u, self.spins, self.sites, siteSpaceTransformation)
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)
        Hamiltonian.__init__(self, [self.spins, self.sites], hubbardMatrix, self.verbose)
        self.sortNsSubspace()

    def sortNsSubspace(self):
        self.sortNs = list()
        nUpEigenvalues = [nsum([1 for digit in self.getOccupationRep(self.fockBasis[fockind])[:int(self.nrOfSingleParticleStates*.5)] if digit == '1']) for fockind in range(self.fockspaceSize)]
        fockstate = 0 # in present basis(most probably N-sorted), eigenvalues evaluated in original basis before
        for blocklength in self.blocksizes:
            found = list()
            sortNsInBlock = list()
            firstBlockEntryNr = fockstate
            for ns in nUpEigenvalues[firstBlockEntryNr:firstBlockEntryNr+blocklength]:
                if ns in found:
                    sortNsInBlock[found.index(ns)].append(fockstate)
                else:
                    sortNsInBlock.append(list())
                    found.append(ns)
                    sortNsInBlock[-1].append(fockstate)
                fockstate += 1
            self.sortNs += sortNsInBlock
        self.blocksizes = [len(block) for block in self.sortNs]
        self.fockBasis = [n for blockn in self.sortNs for n in blockn]
        self.sortNs = coo_matrix(([1]*self.fockspaceSize, (range(self.fockspaceSize), self.fockBasis)), [self.fockspaceSize]*2)
        self.matrix = self.sortNs.dot(self.matrix).dot(self.sortNs.transpose())

    def solve(self):
        Hamiltonian.solve(self)

def setHubbardMatrixWithoutTransf(t, u, spins, orbitals):
    spins = range(len(spins))
    c = AnnihilationOperator([spins, orbitals])
    no = len(orbitals)
    ns = len(spins)
    spininds = range(len(spins))
    ht = [t[i,j] * c[s,i].transpose().dot(c[s,j]) for s in spins for i,j in product(orbitals, orbitals)]
    hu = [.5 * u * c[s1,i].transpose().dot(c[s2, j].transpose()).dot(c[s2, l]).dot(c[s1, k]) for i,j,k,l in product(orbitals, orbitals, orbitals, orbitals) for s1, s2 in product(spins, spins) if s1 != s2]
    return nsum(ht + hu, axis = 0)

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
    ht = [t[i,j] * c[s,i].transpose().dot(c[s,j]) for s in spins for i,j in product(orbitals, orbitals) if t[i,j] != 0]
    hu = [.5 * uMatrix[i, j, k, l, s1, s2] * c[s1,i].transpose().dot(c[s2, j].transpose()).dot(c[s2, l]).dot(c[s1, k]) for i,j,k,l in product(orbitals, orbitals, orbitals, orbitals) for s1, s2 in product(spins, spins) if s1 != s2 and uMatrix[i, j, k, l, s1, s2] != 0]
    return nsum(ht + hu, axis = 0)

def embedV(subspaceVectors, blocksizes):
    assert nsum(blocksizes) == len(subspaceVectors), 'embedding will fail'
    iBlock = 0
    iBlockOrigin = 0
    vectors = zeros([len(subspaceVectors), nsum(blocksizes)])
    for i, v in enumerate(subspaceVectors):
        for j in range(len(v)):
            vectors[i, j + iBlockOrigin] = v[j]
        if i == iBlockOrigin + blocksizes[iBlock] - 1:
            iBlockOrigin += blocksizes[iBlock]
            iBlock += 1
    return vectors
