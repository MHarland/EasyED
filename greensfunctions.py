from numpy import array, sum  as nsum, empty, pi
from time import time
from util import report, scatter_list, allgather_list

class OneParticleGreensFunction(object):
    def __init__(self, verbose = False):
        self.lehmannNominators = dict()
        self.lehmannDenominators = dict()
        self.mesh = None
        self.retardedData = dict()
        self.matsubaraMesh = dict()
        self.matsubaraData = dict()
        self.spectralData = dict()
        self.partitionFunction = None
        self.verbose = verbose

    def setMesh(self, nOmega, omegaMin, omegaMax):
        assert omegaMin <= 0 and omegaMax > 0, 'Choose omega range with 0 in it.'
        self.mesh = array([omegaMin + w*(omegaMax-omegaMin)/float(nOmega-1) for w in range(nOmega)])

    def getMesh(self):
        return self.mesh

    def setRetarded(self, imaginaryOffset = 0.01):
        assert self.partitionFunction != None and self.lehmannNominators != None and self.lehmannDenominators != None, 'Partition Function and Lehmann terms have to be set in advance.'
        report('Calculating one-particle Green\'s function(retarded)...', self.verbose)
        t0 = time()
        for statePair, nominators, denominators in zip(self.lehmannNominators.keys(), self.lehmannNominators.values(), self.lehmannDenominators.values()):
            data_p = list()
            for w in scatter_list(self.mesh):
                terms = [nom/(w + denom + sign * complex(0, imaginaryOffset)) for noms, denom in zip(nominators, denominators) for nom, sign in zip(noms, [+1, +1])]
                data_p.append(nsum(terms/self.partitionFunction))
            self.retardedData.update({statePair: array(allgather_list(data_p))})
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

    def getRetarded(self, statePair):
        if statePair in self.retardedData.keys():
            return self.retardedData[statePair]
        else:
            self.setRetarded()
            return self.retardedData[statePair]

    def setMatsubaraMesh(self, nOmega, beta):
        self.matsubaraMesh = array([(2*n+1)*pi/beta for n in range(nOmega)])

    def getMatsubaraMesh(self):
        return self.matsubaraMesh

    def setMatsubara(self):
        assert self.partitionFunction != None and self.lehmannNominators != None and self.lehmannDenominators != None, 'Partition Function and Lehmann terms have to be set in advance.'
        report('Calculating one-particle Green\'s function(Matsubara)...', self.verbose)
        t0 = time()
        for statePair, nominators, denominators in zip(self.lehmannNominators.keys(), self.lehmannNominators.values(), self.lehmannDenominators.values()):
            data_p = list()
            for wn in scatter_list(self.matsubaraMesh):
                terms = [nom/(denom + complex(0, wn)) for noms, denom in zip(nominators, denominators) for nom in noms]
                data_p.append(nsum(terms/self.partitionFunction))
            self.matsubaraData.update({statePair: array(allgather_list(data_p))})
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
