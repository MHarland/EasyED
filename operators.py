from itertools import product
from numpy import matrix, sum as nsum, binary_repr
from util import diracDelta

class SingleParticleBasis(object):
    def __init__(self, singleParticleBasis = list()):
        self.singleParticleBasis = singleParticleBasis
        self.orderedSingleParticleStates = [a for a in product(*singleParticleBasis)]
        self.nrOfSingleParticleStates = len(self.orderedSingleParticleStates)
        self.fockspaceSize = 2**self.nrOfSingleParticleStates

class AnnihilationOperator(SingleParticleBasis):
    """Fermionic"""
    def __getitem__(self, spState):
        """The single particle state is given by a tuple of quantum numbers."""
        return matrix([[self.element(spState, i, j) for j in range(self.fockspaceSize)] for i in range(self.fockspaceSize)])

    def spsNr(self, *spState):
        """Single particle state number."""
        for i, qnt in enumerate(self.orderedSingleParticleStates):
            if qnt == spState:
                return i
        assert False, 'Single particle state number to state '+str(spState)+' not found.'

    def fsNr(self, *occupationOfSingleParticleStates):
        """Fock state number."""
        return nsum([occ*2**i for i, occ in enumerate(occupationOfSingleParticleStates)])

    def element(self, spState, i, j):
        jbin = binary_repr(j, width = len(self.orderedSingleParticleStates))[::-1]
        spsNr = self.spsNr(*spState)
        nrOfCommutations = nsum([diracDelta('1', jbin[k]) for k in range(spsNr)])
        return (-1)**nrOfCommutations * int(jbin[spsNr]) * diracDelta(i, j - 2**spsNr)

class SuperpositionState(SingleParticleBasis):
    def __init__(self, coefficients, spbasis):
        SingleParticleBasis.__init__(self, spbasis)
        self.coefficients = coefficients

    def getStateAlgebraically(self, thres = 10**(-12)):
        """Fockspace superposition vector."""
        statestr = str()
        for i, coeff in enumerate(self.coefficients):
            if abs(coeff) > thres:
                if len(statestr) > 0:
                    statestr += ' + '
                statestr += str(coeff)
                spstates = list()
                ibin = binary_repr(i, width = self.nrOfSingleParticleStates)
                for j, digit in enumerate(ibin):
                    if digit == '1':
                        statestr += ' c^'+str(self.orderedSingleParticleStates[j])
        return statestr

    def getState(self):
        """Fockspace superposition vector."""
        c = AnnihilationOperator(self.singleParticleBasis)
        matrix = None
        for i, coeff in enumerate(self.coefficients):
            ibin = binary_repr(i, width = self.nrOfSingleParticleStates)
            matrixToAdd = None
            for j, digit in enumerate(ibin):
                if digit == '1' and not matrixToAdd:
                    matrixToAdd = c(*self.orderedSingleParticleStates[j]).H
                elif digit == '1':
                    matrixToAdd = dot(matrixToAdd, c(*self.orderedSingleParticleStates[j]).H)
            if not matrix:
                matrix = coeff * matrixToAdd.copy()
            else:
                matrix += coeff * matrixToAdd
        return matrix
