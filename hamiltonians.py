from itertools import product, izip
from numpy import array, sum as nsum, zeros, argmin, identity, asmatrix, arange, diag, min as nmin
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from time import time

from blocks import BlockMatrix
from operators import SingleParticleBasis, AnnihilationOperator, SuperpositionState
from util import scatter_list, allgather_list, report, contains, getIndex

class Hamiltonian(SingleParticleBasis):
    def __init__(self, singleParticleBasis, matrix, verbose = False):
        self.matrix = matrix
        SingleParticleBasis.__init__(self, singleParticleBasis)
        self.eigenEnergies = None
        self.eigenStates = None # column vectors are eigenstates
        self.blocksizes = [self.fockspaceSize]
        self.verbose = verbose
        self.transformation = coo_matrix(([1]*self.fockspaceSize,
                                          (range(self.fockspaceSize),
                                           range(self.fockspaceSize))))

    def getNEigenvalues(self):
        fockstates = arange(self.fockspaceSize)
        fockstates = self.transformation.dot(fockstates)
        nEigenvalues = [nsum([1 for digit in self.getOccupationRep(fockstate) if digit == '1']) for fockstate in fockstates]
        return nEigenvalues

    def sortSubspaceByPermutation(self, eigenvalues):
        sortedFockstates = list()

        fockstate = 0
        for blocklength in self.blocksizes:
            found = list()
            sortedFockstatesInBlock = list()
            firstBlockState = fockstate
            for eigenvalue in eigenvalues[firstBlockState:firstBlockState+blocklength]:
                if eigenvalue in found:
                    sortedFockstatesInBlock[found.index(eigenvalue)].append(fockstate)
                else:
                    sortedFockstatesInBlock.append(list())
                    found.append(eigenvalue)
                    sortedFockstatesInBlock[-1].append(fockstate)
                fockstate += 1
            sortedFockstates += sortedFockstatesInBlock
        transformationMatrix = coo_matrix(([1]*self.fockspaceSize, (range(self.fockspaceSize), [outState for block in sortedFockstates for outState in block])), [self.fockspaceSize]*2)
        self.blocksizes = [len(block) for block in sortedFockstates]
        self.transformation = transformationMatrix.dot(self.transformation)

    def gatherPermutations(self):
        self.sortSubspaceByPermutation(self.getNEigenvalues())
        
    def solve(self):
        report('Solving the Hamiltonian...', self.verbose)
        t0 = time()
        self.eigenEnergies = list()
        self.eigenStates = list()
        self.gatherPermutations()
        self.matrix = self.transformation.dot(self.matrix).dot(self.transformation.H)
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
        self.eigenEnergies = shiftEnergies(diag(self.transformation.H.dot(self.transformation.transpose().dot(diag(allgather_list(e_scat)).transpose()).transpose())).copy()) # scipy hack for A = U^.D.U
        v = allgather_list(v_scat)
        #self.eigenStates = self.transformation.H.dot(self.transformation.transpose().dot(embedV(v, self.blocksizes)).transpose()) # scipy hack for A = U^.D.U
        self.eigenStates = self.transformation.H.dot(embedV(v, self.blocksizes)).dot(self.transformation)
        self.matrix = self.transformation.H.dot(self.matrix).dot(self.transformation)
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

    def getGroundStateEnergy(self):
        return min(self.eigenEnergies)

    def getGroundSuperpositionState(self, energyResolution = .0001):
        inds = list()
        gss = list()
        for i, e in enumerate(self.eigenEnergies):
            if abs(e - self.getGroundStateEnergy()) < energyResolution:
                inds.append(i)
        for i in inds:
            psi0 = SuperpositionState(self.eigenStates.toarray()[:, i], self.singleParticleBasis)
            gss.append(psi0)
        return gss

    def getGroundStatesAlgebraically(self, energyResolution = .0001):
        groundStates = list()
        for psi in self.getGroundSuperpositionState(energyResolution):
            groundStates.append(psi.getStateAlgebraically())
        return groundStates

    def getGroundStates(self, energyResolution = .0001):
        groundStates = list()
        for psi in self.getGroundSuperpositionState(energyResolution):
            groundStates.append(psi.getStateAlgebraically())
        return groundStates

    def getSpectrum(self, energyResolution = .0001):
        degeneracies = list()
        energies = list()
        for e in self.eigenEnergies:
            if contains(energies, e, energyResolution):
                degeneracies[getIndex(energies, e, energyResolution)] += 1
            else:
                degeneracies.append(1)
                energies.append(e)
        return energies, degeneracies

class Hubbard(Hamiltonian):
    def __init__(self, t, u, siteSpaceTransformation = None, verbose = False):
        self.verbose = verbose
        self.t = array(t)
        self.u = u
        self.spins = ['up', 'dn']
        self.sites = range(len(t))
        self.siteSpaceTransformation = siteSpaceTransformation
        report('Setting up the Hubbard Hamiltonian...', self.verbose)
        t0 = time()
        hubbardMatrix = setHubbardMatrix(self.t, self.u, self.spins, self.sites)
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)
        Hamiltonian.__init__(self, [self.spins, self.sites], hubbardMatrix, self.verbose)

    def getNsEigenvalues(self):
        fockstates = arange(self.fockspaceSize)
        fockstates = self.transformation.dot(fockstates)
        nsEigenvalues = [nsum([1 for digit in self.getOccupationRep(fockstate)[:int(self.nrOfSingleParticleStates*.5)] if digit == '1']) for fockstate in fockstates]
        return nsEigenvalues

    def gatherPermutations(self):
        Hamiltonian.gatherPermutations(self)
        self.sortSubspaceByPermutation(self.getNsEigenvalues())

    def solve(self):
        Hamiltonian.solve(self)

def setHubbardMatrix(t, u, spins, orbitals):
    spins = range(len(spins))
    c = AnnihilationOperator([spins, orbitals])
    no = len(orbitals)
    ns = len(spins)
    spininds = range(len(spins))
    ht = [t[i,j] * c[s,i].H.dot(c[s,j]) for s in spins for i,j in product(orbitals, orbitals) if t[i,j] != 0]
    hu = [.5 * u * c[s1,i].H.dot(c[s1, i]).dot(c[s2, i].H).dot(c[s2, i]) for i in orbitals for s1, s2 in product(spins, spins) if s1 != s2]
    return nsum(ht + hu, axis = 0)
"""
def setHubbardMatrixTransf(t, u, spins, orbitals, siteSpaceTransformation = None): # TODO rm siteSTrafo
    spins = range(len(spins))
    c = AnnihilationOperator([spins, orbitals])
    no = len(orbitals)
    ns = len(spins)
    spininds = range(len(spins))
    uMatrix = zeros([no,no,no,no,ns,ns])
    for i, j, k, l, s1, s2 in product(orbitals,orbitals,orbitals,orbitals,spins,spins):
        if i == k and j == l and i == j and s1 != s2:
            uMatrix[i, j, k, l, s1, s2] = u * .5
    if siteSpaceTransformation != None:
        p = array(siteSpaceTransformation)
        t = p.transpose().dot(t).dot(p)
        temp = uMatrix.copy()
        for i, j, k, l, s1, s2 in product(orbitals,orbitals,orbitals,orbitals,spins,spins):
            uMatrix[i,j,k,l,s1,s2] = nsum([p[i,m] * p[j,n] * temp[m,n,o,q,s1,s2] * p.transpose()[o,l] * p.transpose()[q,k] for m,n,o,q in product(orbitals,orbitals,orbitals,orbitals)], axis = 0)
    ht = [t[i,j] * c[s,i].H.dot(c[s,j]) for s in spins for i,j in product(orbitals, orbitals) if t[i,j] != 0]
    hu = [uMatrix[i, j, k, l, s1, s2] * c[s1,i].H.dot(c[s2, j].H).dot(c[s2, l]).dot(c[s1, k]) for i,j,k,l in product(orbitals, orbitals, orbitals, orbitals) for s1, s2 in product(spins, spins) if s1 != s2 and uMatrix[i, j, k, l, s1, s2] != 0]
    return nsum(ht + hu, axis = 0)
"""
def embedV(subspaceVectors, blocksizes):
    fockspaceSize = nsum(blocksizes)
    assert fockspaceSize == len(subspaceVectors), 'embedding will fail'
    iBlock = 0
    iBlockOrigin = 0
    vectors = zeros([len(subspaceVectors), fockspaceSize])
    x = list()
    y = list()
    data = list()
    for i, v in enumerate(subspaceVectors):
        for j, vj in enumerate(v):
            if vj != 0:
                y.append(j + iBlockOrigin) # TODO understand row/col exchange
                x.append(i)
                data.append(vj)
        if i == iBlockOrigin + blocksizes[iBlock] - 1:
            iBlockOrigin += blocksizes[iBlock]
            iBlock += 1
    return coo_matrix((data, (x,y)), [fockspaceSize]*2)

def shiftEnergies(energies):
    emin = nmin(energies)
    if emin < 0:
        for i in range(len(energies)):
            energies[i] = energies[i] - emin
    return energies
