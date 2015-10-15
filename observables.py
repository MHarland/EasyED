from numpy import array, sum  as nsum, empty, pi
from time import time
from util import report, scatter_list, allgather_list

class StaticObservable(object):
    def __init__(self, operators = dict(), verbose = False):
        self.operators = operators
        self.expectationValue = dict()

    def addOperator(self, state, operator):
        self.operators.update({state: operator})

    def getExpectationValue(self, state):
        assert state in self.operators.keys(), str(state)+' not defined.'
        assert state in self.expectationValue.keys(), str(state)+'\'s expectation value not set.'
        return self.expectationValue[state]

class DynamicObservable(object):
    def __init__(self, operatorPairs = dict(), species = 'fermionic', verbose = False):
        assert species in ['fermionic', 'bosonic'], 'type '+species+' no known.'
        self.operatorPairs = operatorPairs
        self.lehmannNominators = dict()
        self.lehmannDenominators = dict()
        self.zeroFrequencyTerms = dict()
        self.mesh = None
        self.retardedData = dict()
        self.causalData = dict()
        self.matsubaraMesh = None
        self.matsubaraData = dict()
        self.spectralData = dict()
        self.partitionFunction = None
        self.verbose = verbose
        self.customData = dict()
        self.species = species

    def addOperatorPair(self, statePair, operatorPair):
        self.operatorPairs.update({statePair: operatorPair})

    def setMesh(self, nOmega, omegaMin, omegaMax):
        assert omegaMin <= 0 and omegaMax > 0, 'Choose omegaMin <= 0 < omegaMax.'
        self.mesh = array([omegaMin + w*(omegaMax-omegaMin)/float(nOmega-1) for w in range(nOmega)])

    def getMesh(self):
        return self.mesh

    def setRetarded(self, imaginaryOffset):
        assert self.partitionFunction != None and self.lehmannNominators != None and self.lehmannDenominators != None, 'Partition Function and Lehmann terms have to be set in advance.'
        report('Calculating one-particle Green\'s function(retarded)...', self.verbose)
        t0 = time()
        if self.species == 'fermionic':
            self.retardedData.update(lehmannSumDynamic(self.lehmannNominators, self.lehmannDenominators, self.partitionFunction, self.mesh, self.zeroFrequencyTerms, imaginaryOffset, [+1, +1]))
        elif self.species == 'bosonic':
            self.retardedData.update(lehmannSumDynamic(self.lehmannNominators, self.lehmannDenominators, self.partitionFunction, self.mesh, self.zeroFrequencyTerms, imaginaryOffset, [+1, -1]))
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

    def getRetarded(self, statePair):
        if statePair in self.retardedData.keys():
            return self.retardedData[statePair]
        else:
            self.setRetarded()
            return self.retardedData[statePair]

    def setCustom(self, imaginaryOffset, coefficients):
        assert self.partitionFunction != None and self.lehmannNominators != None and self.lehmannDenominators != None, 'Partition Function and Lehmann terms have to be set in advance.'
        report('Calculating custom function...', self.verbose)
        t0 = time()
        self.customData.update(lehmannSumDynamic(self.lehmannNominators, self.lehmannDenominators, self.partitionFunction, self.mesh, self.zeroFrequencyTerms, imaginaryOffset, coefficients))
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

    def getCustom(self, statePair, imaginaryOffset = 0.01, coefficients = [1, 1]):
        if statePair in self.customData.keys():
            return self.customData[statePair]
        else:
            self.setCustom(imaginaryOffset, coefficients)
            return self.customData[statePair]

    def setMatsubaraMesh(self, nOmega, beta):
        if self.species == 'fermionic':
            self.matsubaraMesh = array([complex(0, 2*n+1)*pi/beta for n in range(nOmega)])
        elif self.species == 'bosonic':
            self.matsubaraMesh = array([complex(0, 2*n)*pi/beta for n in range(nOmega)])

    def getMatsubaraMesh(self):
        return self.matsubaraMesh

    def setMatsubara(self):
        assert self.partitionFunction != None and self.lehmannNominators != None and self.lehmannDenominators != None, 'Partition Function and Lehmann terms have to be set in advance.'
        report('Calculating one-particle Green\'s function(Matsubara)...', self.verbose)
        t0 = time()
        if self.species == 'fermionic':
            self.matsubaraData.update(lehmannSumDynamic(self.lehmannNominators, self.lehmannDenominators, self.partitionFunction, self.matsubaraMesh, self.zeroFrequencyTerms, 0, [+1, +1]))
        if self.species == 'bosonic':
            self.matsubaraData.update(lehmannSumDynamic(self.lehmannNominators, self.lehmannDenominators, self.partitionFunction, self.matsubaraMesh, self.zeroFrequencyTerms, 0, [+1, -1]))
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

    def getMatsubara(self, statePair):
        if statePair in self.matsubaraData.keys():
            return self.matsubaraData[statePair]
        else:
            self.setMatsubara()
            return self.matsubaraData[statePair]

    def getSpectralFunction(self, statePair):
        """use explicit setRetarded to change the imaginary Offset"""
        return array([-1/pi *gret.imag for gret in self.getRetarded(statePair)])
    
def lehmannSumDynamic(elementProducts, energyDifferences, partitionFunction, mesh, zeroFrequencyTerms, imaginaryOffset, nominatorCoefficients):
    if type(imaginaryOffset) != list:
        imaginaryOffset = [imaginaryOffset]*2
    result = dict()
    for statePair, nominators, denominators in zip(elementProducts.keys(), elementProducts.values(), energyDifferences.values()):
        data_p = list()
        for w in scatter_list(mesh):
            terms = [coeff * nom/(w + denom + complex(0, offset)) for noms, denom in zip(nominators, denominators) for nom, coeff, offset in zip(noms, nominatorCoefficients, imaginaryOffset)]
            if len(zeroFrequencyTerms.values()) > 0 and w == 0:
                terms += zeroFrequencyTerms[statePair]
            data_p.append(nsum(terms/partitionFunction))
        result.update({statePair: array(allgather_list(data_p))})
    return result
