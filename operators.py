from itertools import product
from numpy import matrix, sum as nsum, binary_repr, dot, zeros, identity, sqrt, bitwise_xor

from blocks import BlockMatrix
from util import diracDelta

class SingleParticleBasis(object):
    def __init__(self, singleParticleBasis = list()):
        self.singleParticleBasis = singleParticleBasis
        self.orderedSingleParticleStates = [a for a in product(*singleParticleBasis)]
        self.nrOfSingleParticleStates = len(self.orderedSingleParticleStates)
        self.fockspaceSize = 2**self.nrOfSingleParticleStates

    def getState(self, fockspaceNr):
        statestr = str()
        binCode = binary_repr(fockspaceNr, width = self.nrOfSingleParticleStates)
        for i, digit in enumerate(binCode):
            if digit == '1':
                statestr += ' c^'+str(self.orderedSingleParticleStates[i])
        return statestr

    def getSingleParticleBasis(self):
        return self.singleParticleBasis

    def getSingleParticleStateNr(self, *spState):
        """Single particle state number."""
        return self.orderedSingleParticleStates.index(spState)

    def getFockspaceNr(self, *occupationOfSingleParticleStates):
        """Fock state number."""
        return nsum([occ*2**i for i, occ in enumerate(occupationOfSingleParticleStates)])

    def getBinaryRep(self, fockspaceNr):
        return binary_repr(fockspaceNr, width = len(self.orderedSingleParticleStates))[::-1]

class AnnihilationOperator(SingleParticleBasis):
    """Fermionic"""
    def __init__(self, singleParticleBasis, sortNSubspaces = True):
        SingleParticleBasis.__init__(self, singleParticleBasis)
        if sortNSubspaces:
            self.sortN = True
            self.blocksizes = list()
            self.sortByN = list()
            nEigenvalues = [nsum([1 for digit in self.getBinaryRep(fockind) if digit == '1']) for fockind in range(self.fockspaceSize)]
            found = list()
            for fockspaceNr, n in enumerate(nEigenvalues):
                if n in found:
                    self.sortByN[found.index(n)].append(fockspaceNr)
                else:
                    self.sortByN.append(list())
                    found.append(n)
                    self.sortByN[-1].append(fockspaceNr)
            self.blocksizes = [len(block) for block in self.sortByN]
            self.sortByN = [n for blockn in self.sortByN for n in blockn]

    def __getitem__(self, spState):
        """The single particle state is given by a tuple of quantum numbers. If sortNSubspaces, then it returns a list of blocks, i.e. a list of 2D arrays"""
        spsnr = self.getSingleParticleStateNr(*spState)
        targetStates = list(bitwise_xor([range(self.fockspaceSize)], 2**spsnr)[0]) # array better
        signs = [(-1)**nsum([diracDelta('1', self.getBinaryRep(k)) for k in range(spsnr)])]# replace diracD by if after range
        print 'ttt', signs
        targetStatesT = self.nBlockTransformation(targetStates)
        signsT = self.nBlockTransformation(signs)
        instatesT = self.nBlockTransformation(range(self.fockspaceSize))
        matrix = BlockMatrix(self.blocksizes)
        for j, targetState, sign in zip(instatesT, targetStatesT, signsT):
            matrix.setFullMatrixEntry(targetState, j, sign)
        return matrix

    def nBlockTransformation(self, states):
        return [states[j] for j in self.sortByN]

    def nBlockBackTransformation(self, states):
        return [states[states.index(j)] for j in self.sortByN]

    def element(self, spState, i, j):
        jbin = binary_repr(j, width = len(self.orderedSingleParticleStates))[::-1]
        spsNr = self.spsNr(*spState)
        nrOfCommutations = nsum([diracDelta('1', jbin[k]) for k in range(spsNr)])
        return (-1)**nrOfCommutations * int(jbin[spsNr]) * diracDelta(i, j - 2**spsNr)

class SuperpositionState(SingleParticleBasis):
    def __init__(self, coefficients, spbasis):
        SingleParticleBasis.__init__(self, spbasis)
        self.coefficients = coefficients

    def getState(self, thres = 10**(-12)):
        """Fockspace superposition vector."""
        statestr = str()
        for i, coeff in enumerate(self.coefficients):
            if abs(coeff) > thres:
                if len(statestr) > 0:
                    statestr += ' + '
                statestr += str(coeff)
                spstates = list()
                statestr += self.getState(i)
        return statestr

    def getMatrix(self):
        """Fockspace superposition vector."""
        c = AnnihilationOperator(self.singleParticleBasis)
        stateMatrix = zeros([c.fockspaceSize, c.fockspaceSize])
        for i, coeff in enumerate(self.coefficients):
            ibin = binary_repr(i, width = self.nrOfSingleParticleStates)
            matrixToAdd = identity(c.fockspaceSize)
            for j, digit in enumerate(ibin):
                if digit == '1':
                    matrixToAdd = dot(matrixToAdd, c[self.orderedSingleParticleStates[j]].H)
            stateMatrix += coeff * matrixToAdd
        return matrix(stateMatrix)
