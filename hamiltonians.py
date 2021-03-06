from itertools import product, izip
from numpy import array, sum as nsum, zeros, argmin, identity, asmatrix, arange, diag, min as nmin, argsort, sort
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from time import time

from blocks import BlockMatrix
from operators import SingleParticleBasis, AnnihilationOperator, SuperpositionState
from util import scatter_list, allgather_list, report, contains, getIndex, equals

class Hamiltonian(SingleParticleBasis):
    def __init__(self, singleParticleBasis, matrix, verbose = False, all_real = True):
        self.matrix = matrix
        SingleParticleBasis.__init__(self, singleParticleBasis)
        self.eigenEnergies = None
        self.eigenStates = None # column vectors are eigenstates
        self.blocksizes = [self.fockspaceSize]
        self.verbose = verbose
        self.transformation = coo_matrix(([1]*self.fockspaceSize,
                                          (range(self.fockspaceSize),
                                           range(self.fockspaceSize))))
        self.energyShift = None
        self.all_real = all_real
        self.matrix_sorted = None

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
        
    def solve(self, store_sorted_matrix = False):
        report('Solving the Hamiltonian...', self.verbose)
        t0 = time()
        self.eigenEnergies = list()
        self.eigenStates = list()
        self.gatherPermutations()
        self.matrix = self.transformation.dot(self.matrix).dot(self.transformation.H)
        if store_sorted_matrix:
            self.matrix_sorted = self.matrix.copy()
        hBlocks = BlockMatrix(self.blocksizes, self.all_real)
        rows, cols = self.matrix.nonzero()
        for val, i, j in izip(self.matrix.data, rows, cols):
            hBlocks[i, j] = val
        #print dir(hBlocks)
        #print hBlocks.getBlock(4)
        hBlocks_scat = scatter_list(hBlocks.datablocks)
        e_scat = list()
        v_scat = list()
        for block in hBlocks_scat:
            e, v = eigh(block)
            e_scat += list(e)
            v_scat += list(v)
        self.eigenEnergies, self.energyShift = shiftEnergies(diag(self.transformation.H.dot(self.transformation.transpose().dot(diag(allgather_list(e_scat)).transpose()).transpose())).copy()) # scipy hack for A = U^.D.U
        #self.eigenEnergies = diag(self.transformation.H.dot(self.transformation.transpose().dot(diag(allgather_list(e_scat)).transpose()).transpose())).copy() # scipy hack for A = U^.D.U
        v = allgather_list(v_scat)
        self.eigenStates = self.transformation.H.dot(embedV(v, self.blocksizes)).dot(self.transformation)
        self.matrix = self.transformation.H.dot(self.matrix).dot(self.transformation)
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

    def getGroundStateEnergy(self):
        #return min(self.eigenEnergies)
        return self.energyShift

    def getGroundSuperpositionState(self, energyResolution = 10**(-10)):
        inds = list()
        gss = list()
        for i, e in enumerate(self.eigenEnergies):
            if abs(e - 0) < energyResolution:
                inds.append(i)
        for i in inds:
            psi0 = SuperpositionState(self.eigenStates.toarray()[:, i], self.singleParticleBasis)
            gss.append(psi0)
        return gss

    def getSuperpositionStatesEnergySorted(self, energyResolution = 10**(-10), threshold = .0001):
        states = list()
        states.append(list())
        inds = argsort(self.eigenEnergies)
        for j, i in enumerate(inds):
            e = self.eigenEnergies[i]
            if j > 0:
                if not abs(e-self.eigenEnergies[inds[j-1]]) < energyResolution:
                    states.append(list())
            psi0 = SuperpositionState(self.eigenStates.toarray()[:, i], self.singleParticleBasis)
            states[-1].append(psi0)
        return states

    def getStatesEnergySortedAlgebraically(self, energyResolution = 10**(-10), threshold = .0001, n_coeff = None):
        states = list()
        for energyGroup in self.getSuperpositionStatesEnergySorted(energyResolution):
            states.append(list())
            for psi in energyGroup:
                if n_coeff == None:
                    states[-1].append(psi.getStateAlgebraically(threshold))
                else:
                    states[-1].append(psi.getStateAlgebraicallyByCoefficient(n_coeff = n_coeff))
        return states

    def getGroundStatesAlgebraically(self, energyResolution = 10**(-10), threshold = .0001):
        groundStates = list()
        for psi in self.getGroundSuperpositionState(energyResolution):
            groundStates.append(psi.getStateAlgebraically(threshold))
        return groundStates

    def getGroundStates(self, energyResolution = 10**(-10), threshold = .0001):
        groundStates = list()
        for psi in self.getGroundSuperpositionState(energyResolution):
            groundStates.append(psi)
        return groundStates

    def getEnergies(self, energyResolution = 10**(-10)):
        energies = list()
        for i, e in enumerate(sort(self.eigenEnergies)):
            if i > 0:
                if abs(e-eold) > energyResolution:
                    energies.append(e)
            else:
                energies.append(e)
            eold = e
        return energies

    def getSpectrum(self, energyResolution = 10**(-10)):
        degeneracies = list()
        energies = list()
        for e in self.eigenEnergies:
            if contains(energies, e, energyResolution):
                degeneracies[getIndex(energies, e, energyResolution)] += 1
            else:
                degeneracies.append(1)
                energies.append(e)
        return energies, degeneracies

    def getSpectrumEnergySorted(self, energyResolution = 10**(-10)):
        energies, degeneracies = self.getSpectrum(energyResolution)
        inds = argsort(energies)
        energiesSorted = []
        degeneraciesSorted = []
        for ind in inds:
            energiesSorted.append(energies[ind])
            degeneraciesSorted.append(degeneracies[ind])
        return energiesSorted, degeneraciesSorted


class Hubbard(Hamiltonian):
    """
    sp_basis must have spins first, then sites
    """
    def __init__(self, t, u, siteSpaceTransformation = None, transformationLabels = None, verbose = False, n_is_quantumnr = True, sz_is_quantumnr = True, parity_is_quantumnr = False, sp_basis = None):
        self.verbose = verbose
        self.t = array(t)
        self.u = u
        self.spins = ['up', 'dn']
        self.sites = range(len(t))
        self.siteSpaceTransformation = siteSpaceTransformation
        self.n_is_quantumnr = n_is_quantumnr
        self.sz_is_quantumnr = sz_is_quantumnr
        self.parity_is_quantumnr = parity_is_quantumnr
        report('Setting up the Hubbard Hamiltonian...', self.verbose)
        t0 = time()
        hubbardMatrix = setHubbardMatrix(self.t, self.u, self.spins, self.sites, siteSpaceTransformation)
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)
        if transformationLabels == None:
            orbitals = self.sites
        else:
            orbitals =  transformationLabels
        if sp_basis is None:
            sp_basis = [self.spins, orbitals]
        Hamiltonian.__init__(self, sp_basis, hubbardMatrix, self.verbose)

    def getSzEigenvalues(self):
        fockstates = arange(self.fockspaceSize)
        fockstates = self.transformation.dot(fockstates)
        szEigenvalues = [.5*(nsum([1 for digit in self.getOccupationRep(fockstate)[:int(self.nrOfSingleParticleStates*.5)] if digit == '1']) - nsum([1 for digit in self.getOccupationRep(fockstate)[int(self.nrOfSingleParticleStates*.5):] if digit == '1'])) for fockstate in fockstates]
        return szEigenvalues

    def getParityEigenvalues(self):
        fockstates = arange(self.fockspaceSize)
        fockstates = self.transformation.dot(fockstates)
        pEigenvalues = [nsum([1 for digit in self.getOccupationRep(fockstate) if digit == '1'])%2 for fockstate in fockstates]
        return pEigenvalues
    
    def gatherPermutations(self):
        if self.n_is_quantumnr:
            Hamiltonian.gatherPermutations(self)
        if self.parity_is_quantumnr:
            self.sortSubspaceByPermutation(self.getParityEigenvalues())
        if self.sz_is_quantumnr:
            self.sortSubspaceByPermutation(self.getSzEigenvalues())

"""
    def solve(self, *args, **kwargs):
        Hamiltonian.solve(self, *args, **kwargs)
"""

"""
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

class Hubbard2(Hubbard):
    """
    This class uses the single particle basis sites x spins, which is uncommon, but nessecary for the calculation of the reduced density matrix
    sp_basis must have sites first then spins
    """
    def __init__(self, t, u, verbose = False, n_is_quantumnr = True, sz_is_quantumnr = True, parity_is_quantumnr = False, sp_basis = None):
        self.verbose = verbose
        self.t = array(t)
        self.u = u
        self.spins = range(2)#['up', 'dn']
        self.sites = range(len(t))
        self.n_is_quantumnr = n_is_quantumnr
        self.sz_is_quantumnr = sz_is_quantumnr
        self.parity_is_quantumnr = parity_is_quantumnr
        no = len(self.sites)
        ns = len(self.spins)
        if sp_basis is None:
            c = AnnihilationOperator([self.sites, self.spins])
            sp_basis = [self.sites, self.spins]
        else:
            c = AnnihilationOperator([sp_basis[0], sp_basis[1]])
        uMatrix = zeros([no,no,no,no,ns,ns])
        for i, j, k, l, s1, s2 in product(*([self.sites]*4 + [self.spins]*2)):
            if i == k and j == l and i == j and s1 != s2:
                uMatrix[i, j, k, l, s1, s2] = u * .5
        t = array(t)
        ht = [t[i,j] * c[i,s].H.dot(c[j,s]) for s in self.spins for i,j in product(*[self.sites]*2) if t[i,j] != 0]
        hu = [uMatrix[i, j, k, l, s1, s2] * c[i,s1].H.dot(c[j,s2].H).dot(c[l,s2]).dot(c[k,s1]) for i,j,k,l in product(*[self.sites]*4) for s1, s2 in product(*[self.spins]*2) if s1 != s2 and uMatrix[i, j, k, l, s1, s2] != 0]
        hubbardMatrix = nsum(ht + hu, axis = 0)
        Hamiltonian.__init__(self, sp_basis, hubbardMatrix, self.verbose)

    def getSzEigenvalues(self):
        fockstates = arange(self.fockspaceSize)
        fockstates = self.transformation.dot(fockstates)
        szEigenvalues = [.5*(nsum([1 for digit in self.getOccupationRep(fockstate)[::2] if digit == '1']) - nsum([1 for digit in self.getOccupationRep(fockstate)[1::2] if digit == '1'])) for fockstate in fockstates]
        return szEigenvalues

    def gatherPermutations(self):
        if self.sz_is_quantumnr:
            self.sortSubspaceByPermutation(self.getSzEigenvalues())
        if self.n_is_quantumnr:
            Hamiltonian.gatherPermutations(self)
        if self.parity_is_quantumnr:
            self.sortSubspaceByPermutation(self.getParityEigenvalues())
    
def setHubbardMatrix(t, u, spins, orbitals, siteSpaceTransformation = None): # TODO rm siteSTrafo
    spins = range(len(spins))
    c = AnnihilationOperator([spins, orbitals])
    no = len(orbitals)
    ns = len(spins)
    spininds = range(len(spins))
    uMatrix = zeros([no,no,no,no,ns,ns])
    for i, j, k, l, s1, s2 in product(orbitals,orbitals,orbitals,orbitals,spins,spins):
        if i == k and j == l and i == j and s1 != s2:
            uMatrix[i, j, k, l, s1, s2] = u * .5
    if siteSpaceTransformation is not None:
        p = array(siteSpaceTransformation)
        t = p.transpose().dot(t).dot(p)
        temp = uMatrix.copy()
        for i, j, k, l, s1, s2 in product(orbitals,orbitals,orbitals,orbitals,spins,spins):
            uMatrix[i,j,k,l,s1,s2] = nsum([p[i,m] * p[j,n] * temp[m,n,o,q,s1,s2] * p.transpose()[o,l] * p.transpose()[q,k] for m,n,o,q in product(orbitals,orbitals,orbitals,orbitals)], axis = 0)
    ht = [t[i,j] * c[s,i].H.dot(c[s,j]) for s in spins for i,j in product(orbitals, orbitals) if t[i,j] != 0]
    hu = [uMatrix[i, j, k, l, s1, s2] * c[s1,i].H.dot(c[s2, j].H).dot(c[s2, l]).dot(c[s1, k]) for i,j,k,l in product(orbitals, orbitals, orbitals, orbitals) for s1, s2 in product(spins, spins) if s1 != s2 and uMatrix[i, j, k, l, s1, s2] != 0]
    return nsum(ht + hu, axis = 0)

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
            if not equals(vj, 0):
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
    return energies, emin
