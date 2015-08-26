from itertools import product
from numpy import trace, exp, sum as nsum, identity, where, array
from scipy.optimize import bisect
from time import time
from greensfunctions import OneParticleGreensFunction
from operators import AnnihilationOperator
from util import scatter_list, sumScatteredLists, report, allgather_list

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
        self.g1 = OneParticleGreensFunction(self.verbose)

    def calcOccupation(self, singleParticleState = None):
        c = AnnihilationOperator(self.singleParticleBasis)
        e = self.getEnergyEigenvalues()
        ve = self.getEnergyEigenstates().toarray()
        z = self.getPartitionFunction()
        fockstates = range(self.fockspaceSize)
        if singleParticleState == None:
            states = self.orderedSingleParticleStates
        else:
            states = [singleParticleState]

        report('Calculating Occupation...', self.verbose)
        t0 = time()
        for state in states:
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

    def getEnergyEigenstates(self):
        if self.energyEigenstates != None:
            return self.energyEigenstates
        else:
            self.calcEigensystem()
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

    def calcG1(self, singleParticleStatePairs = None):
        c = AnnihilationOperator(self.singleParticleBasis)
        e = self.getEnergyEigenvalues()
        ve = self.getEnergyEigenstates().toarray()
        fockstates = range(self.fockspaceSize)
        if singleParticleStatePairs == None:
            statePairs = [(state1, state2) for state1, state2 in product(self.orderedSingleParticleStates, self.orderedSingleParticleStates)]
        else:
            statePairs = singleParticleStatePairs

        self.g1.partitionFunction = self.getPartitionFunction()

        report('Calculating one-particle Green\'s function...', self.verbose)
        t0 = time()
        for statePair in statePairs:
            c_state = c[statePair[0]].toarray()
            c_state_dag = c[statePair[1]].H.toarray()
            n_inds_p = array(scatter_list(fockstates))
            nominators_p = []
            denominators_p = []
            for n in n_inds_p:
                n_c_state = ve[:,n].dot(c_state)
                c_state_dag_n = c_state_dag.dot(ve[:,n])
                expn = exp(-self.beta*e[n])
                for m in fockstates:
                    el1 = n_c_state.dot(ve[:,m])
                    if el1 != 0:
                        el2 = ve[:,m].dot(c_state_dag_n)
                        if el2 != 0:
                            expm = exp(-self.beta*e[m])
                            el = el1*el2
                            nominators_p.append([])
                            nominators_p[-1].append(expn*el)
                            nominators_p[-1].append(expm*el)
                            denominators_p.append(e[n]-e[m])
            self.g1.lehmannNominators.update({statePair: allgather_list(nominators_p)})
            self.g1.lehmannDenominators.update({statePair: allgather_list(denominators_p)})
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

class GrandcanonicalEnsemble(CanonicalEnsemble):
    """use canonical with mu in H, TODO transform mu"""
    def __init__(self, hamiltonian, beta, mu, verbose = True):
        self.mu = mu
        CanonicalEnsemble.__init__(self, hamiltonian, beta, verbose)
        c = AnnihilationOperator(self.singleParticleBasis)
        muMatrix = mu * nsum([c[orb].H.dot(c[orb]) for orb in self.orderedSingleParticleStates], axis = 0)
        self.hamiltonian.matrix = self.hamiltonian.matrix - muMatrix

    def setMu(self, filling, muMin, muMax, muTol = .001, maxiter = 100):
        c = AnnihilationOperator(self.singleParticleBasis)
        nMatrix = nsum([c[orb].H.dot(c[orb]) for orb in self.orderedSingleParticleStates], axis = 0)
        self.hamiltonian.matrix = self.hamiltonian.matrix + self.mu * nMatrix
        def fillingFunction(muTrial):
            self.hamiltonian.matrix = self.hamiltonian.matrix - muTrial * nMatrix
            self.calcEigensystem()
            self.calcOccupation()
            fillingTrial = self.getTotalOccupation()
            self.hamiltonian.matrix = self.hamiltonian.matrix + muTrial * nMatrix
            report('Filling(mu='+str(muTrial)+') = '+str(fillingTrial), self.verbose)
            return fillingTrial - filling
        mu = bisect(fillingFunction, muMin, muMax, xtol = muTol, maxiter = maxiter)
        self.mu = mu
        self.hamiltonian.matrix = self.hamiltonian.matrix - mu * nMatrix
        report('Chemical potential set to '+str(mu), self.verbose)

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
