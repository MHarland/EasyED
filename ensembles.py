from itertools import product
from numpy import trace, exp, sum as nsum, identity, where, array
from time import time
from operators import AnnihilationOperator
from util import scatter_list, sumScatteredLists, report

class CanonicalEnsemble(object):
    def __init__(self, hamiltonian, beta, verbose = True):
        self.hamiltonian = hamiltonian
        self.singleParticleBasis = hamiltonian.singleParticleBasis
        self.orderedSingleParticleStates = hamiltonian.orderedSingleParticleStates
        self.fockspaceSize = hamiltonian.fockspaceSize
        self.partitionFunction = None
        self.beta = beta
        self.energyEigenstates = None
        self.energyEigenvalues = None
        self.occupation = dict()
        self.verbose = verbose

    def calcOccupation(self):
        c = AnnihilationOperator(self.singleParticleBasis)
        e = self.getEnergyEigenvalues()
        ve = self.getEnergyEigenstates().toarray()
        z = self.getPartitionFunction()
        inds = self.getIndsForSum()
        fockstates = range(self.fockspaceSize)

        report('Calculating Occupation...', self.verbose)
        t0 = time()
        for state in self.orderedSingleParticleStates:
            c_state_dag = c[state].H.toarray()
            n_inds_p = array(scatter_list(fockstates))
            terms_p = list()
            for n in n_inds_p:
                c_state_dag_n = c_state_dag.dot(ve[:,n])
                expn = exp(-self.beta*e[n])
                for m in fockstates:
                    el = ve[:,m].dot(c_state_dag_n)
                    if el != 0:
                        terms_p.append(expn*el*el.conjugate())
            self.occupation.update({state: sumScatteredLists(terms_p)/z})
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

    def getTotalOccupation(self):
        return nsum(self.occupation.values())

    def calcPartitionFunction(self):
        report('Calculating Partition Function...', self.verbose)
        t0 = time()
        self.partitionFunction = nsum([exp(-self.beta*eva) for eva in self.getEnergyEigenvalues()])
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

    def getPartitionFunction(self):
        if self.partitionFunction != None:
            return self.partitionFunction
        else:
            self.calcPartitionFunction()
            return self.partitionFunction

    def getIndsForSum(self):
        inds = range(self.fockspaceSize)
        return [(m, n) for m, n in product(inds, inds)]

    def getEnergyEigenstates(self):
        if self.energyEigenstates != None:
            return self.energyEigenstates
        else:
            self.self.calcEigensystem()
            return self.energyEigenstates

    def getEnergyEigenvalues(self):
        if self.energyEigenvalues != None:
            return self.energyEigenvalues
        else:
            self.calcEigensystem()
            return self.energyEigenvalues

    def calcEigensystem(self):
        self.hamiltonian.solve()
        self.energyEigenstates = self.hamiltonian.eigenStates
        self.energyEigenvalues = self.hamiltonian.eigenEnergies

class GrandcanonicalEnsemble(CanonicalEnsemble):
    """use canonical with mu in H, TODO transform mu"""
    def __init__(self, hamiltonian, beta, mu):
        CanonicalEnsemble.__init__(self, hamiltonian, beta)
        c = AnnihilationOperator(self.singleParticleBasis)
        muMatrix = mu * nsum([c[orb].H.dot(c[orb]) for orb in self.orderedSingleParticleStates], axis = 0)
        self.hamiltonian.matrix = self.hamiltonian.matrix - muMatrix

class MicrocanonicalEnsemble(CanonicalEnsemble):
    def __init__(self, hamiltonian):
        CanonicalEnsemble.__init__(self, hamiltonian, 0)

    #def getIndsForSum(self):
    #    inds, = where(self.energyEigenvalues == self.energyEigenvalues.min())
    #    print inds
    #    return [(m, n) for m, n in product(inds, inds)]

    def calcPartitionFunction(self):
        #inds, = where(self.energyEigenvalues == self.energyEigenvalues.min())
        #self.partitionFunction = nsum([1 for eva in inds])
        self.partitionFunction = self.fockspaceSize
