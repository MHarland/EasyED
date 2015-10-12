from itertools import product
from numpy import matrix, sum as nsum, binary_repr, dot, zeros, identity, sqrt, bitwise_xor, sort, argsort
from scipy.sparse import coo_matrix

class SingleParticleBasis(object):
    def __init__(self, singleParticleBasis = list()):
        self.singleParticleBasis = singleParticleBasis
        self.orderedSingleParticleStates = [a for a in product(*singleParticleBasis)]
        self.nrOfSingleParticleStates = len(self.orderedSingleParticleStates)
        self.fockspaceSize = 2**self.nrOfSingleParticleStates

    def getStateAlgebraically(self, fockspaceNr):
        statestr = str()
        for i, digit in enumerate(self.getOccupationRep(fockspaceNr)):
            if digit == '1':
                if i > 0:
                    statestr += ' '
                statestr += 'c^'+str(self.orderedSingleParticleStates[i])
        return statestr

    def getSingleParticleBasis(self):
        return self.singleParticleBasis

    def getSingleParticleStateNr(self, *spState):
        """Single particle state number."""
        return self.orderedSingleParticleStates.index(spState)

    def getFockspaceNr(self, occupationOfSingleParticleStates = None, singleParticleStateNr = None, singleParticleState = None):
        """Fock state number."""
        if occupationOfSingleParticleStates != None:
            return nsum([int(occ)*2**i for i, occ in enumerate(occupationOfSingleParticleStates)])
        elif singleParticleStateNr != None:
            return 2**singleParticleStateNr
        elif singleParticleState != None:
            return 2**self.getSingleParticleStateNr(*singleParticleState)
        else:
            assert False, 'Need parameter.'

    def getOccupationRep(self, fockstateNr = None, singleParticleStateNr = None, singleParticleState = None):
        if fockstateNr != None:
            return binary_repr(fockstateNr, width = len(self.orderedSingleParticleStates))[::-1]
        elif singleParticleStateNr != None:
            occstr = str()
            for i in range(self.nrOfSingleParticleStates):
                if i == singleParticleStateNr:
                    occstr += '1'
                else:
                    occstr += '0'
            return occstr

        elif singleParticleState != None:
            occstr = str()
            for i in range(self.nrOfSingleParticleStates):
                if i == self.getSingleParticleStateNr(*singleParticleState):
                    occstr += '1'
                else:
                    occstr += '0'
            return occstr

class AnnihilationOperator(SingleParticleBasis):
    """Fermionic"""

    def __getitem__(self, spState):
        """The single particle state is given by a tuple of quantum numbers. Returns scipy.sparse.coo_matrix"""
        instates = list()
        outstates = list()
        spStateOR = self.getOccupationRep(singleParticleState = spState)
        for fockStateNr in range(self.fockspaceSize):
            instateOR = self.getOccupationRep(fockStateNr)
            if instateOR[self.orderedSingleParticleStates.index(spState)] == '1':
                instates.append(fockStateNr)
                outstates.append(self.getFockspaceNr(annihilateOccRep(spStateOR, instateOR)))
        signs = [(-1)**nsum([1 for k in range(self.getSingleParticleStateNr(*spState)) if self.getOccupationRep(fockstateNr)[k] == '1']) for fockstateNr in instates]
        return coo_matrix((signs, (outstates, instates)), [self.fockspaceSize]*2)

class SuperpositionState(SingleParticleBasis):
    def __init__(self, coefficients, spbasis):
        SingleParticleBasis.__init__(self, spbasis)
        self.coefficients = coefficients

    def getStateAlgebraically(self, thres = .0001, linebreak = True):
        statestr = str()
        for i in argsort(abs(self.coefficients))[::-1]:
            coeff = self.coefficients[i]
            if abs(coeff) >= thres:
                if len(statestr) > 0 and not linebreak:
                    statestr += ' '
                if coeff > 0:
                    statestr += '+'
                statestr += str(coeff)[:7]
                spstates = list()
                statestr += SingleParticleBasis.getStateAlgebraically(self, i)
                if linebreak:
                    statestr += '\n'
        return statestr

    def getStateAlgebraicallyByCoefficient(self, energyResolution = 10**(-10), n_coeff = 1):
        threshold = sort(abs(self.coefficients))[-n_coeff] - energyResolution
        return self.getStateAlgebraically(threshold)

    def getQuantumNumber(self, operator):
        return self.coefficients.dot(operator.dot(self.coefficients))

def annihilateOccRep(spsToAnnihilate, fockstate):
    newState = str()
    for i, occ in enumerate(fockstate):
        if spsToAnnihilate[i] == '0':
            newState += occ
        else:
            newState += '0'
    return newState
