from itertools import product
from numpy import trace, exp, sum as nsum, identity, where, array, ndarray
from scipy.optimize import bisect
from time import time
from observables import StaticObservable, DynamicObservable
from operators import AnnihilationOperator
from util import scatter_list, sumScatteredLists, report, allgather_list, equals

class CanonicalEnsemble(object):
    def __init__(self, hamiltonian, beta, verbose = False):
        self.hamiltonian = hamiltonian
        self.singleParticleBasis = hamiltonian.singleParticleBasis
        self.orderedSingleParticleStates = hamiltonian.orderedSingleParticleStates
        self.fockspaceSize = hamiltonian.fockspaceSize
        self.partitionFunction = None
        self.beta = beta
        self.energyEigenstates = None
        self.energyEigenvalues = None
        self.verbose = verbose
        self.occupation = StaticObservable(verbose = self.verbose)
        self.g1 = DynamicObservable(verbose = self.verbose)

    def setLehmannSumStatic(self, staticObservable):
        assert type(staticObservable) == StaticObservable, 'StaticObservable expected.'
        e = self.getEnergyEigenvalues()
        ve = self.getEnergyEigenstates().toarray()
        z = self.getPartitionFunction()
        fockstates = range(self.fockspaceSize)
        for index, matrix in staticObservable.operators.items():
            if type(matrix) != ndarray:
                m = matrix.toarray()
            else:
                m = matrix
            n_inds_p = array(scatter_list(fockstates))
            terms_p = list()
            for n in n_inds_p:
                el = ve[:,n].dot(m.dot(ve[:,n]))
                expn = exp(-self.beta*e[n])
                if not equals(el, 0) and not equals(expn, 0):
                    terms_p.append(expn * el)
            staticObservable.expectationValue.update({index: sumScatteredLists(terms_p)/z})

    def calcOccupation(self, singleParticleState = None):
        c = AnnihilationOperator(self.singleParticleBasis)
        if singleParticleState == None:
            states = self.orderedSingleParticleStates
        else:
            states = [singleParticleState]
        report('Calculating Occupation...', self.verbose)
        t0 = time()
        for state in states:
            n_state = c[state].H.dot(c[state])
            self.occupation.addOperator(state, n_state)
        self.setLehmannSumStatic(self.occupation)
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

    def getTotalOccupation(self):
        return nsum(self.occupation.expectationValue.values())

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

    def setLehmannTermsDynamic(self, dynamicObservable):
        assert type(dynamicObservable) == DynamicObservable, 'Type DynamicObservable expected.'
        if dynamicObservable.species == 'fermionic':
            self.setLehmannTermsDynamicFermionic(dynamicObservable)
        elif dynamicObservable.species == 'bosonic':
            self.setLehmannTermsDynamicBosonic(dynamicObservable)

    def setLehmannTermsDynamicFermionic(self, dynamicObservable):
        for statePair, operatorPair in dynamicObservable.operatorPairs.items():
            operator1 = operatorPair[0]
            operator2 = operatorPair[1]
            if type(operator1) != ndarray:
                operator1 = operator1.toarray()
                operator2 = operator2.toarray()
            e = self.getEnergyEigenvalues()
            ve = self.getEnergyEigenstates().toarray()
            fockstates = range(self.fockspaceSize)

            n_inds_p = array(scatter_list(fockstates))
            nominators_p = []
            denominators_p = []
            for n in n_inds_p:
                n_op1 = ve[:,n].dot(operator1)
                op2_n = operator2.dot(ve[:,n])
                expn = exp(-self.beta*e[n])
                for m in fockstates:
                    if dynamicObservable.species == 'bosonic':
                        if equals(e[n], e[m]):
                            continue
                    el1 = n_op1.dot(ve[:,m])
                    if not equals(el1, 0):
                        el2 = ve[:,m].dot(op2_n)
                        if not equals(el2, 0):
                            expm = exp(-self.beta*e[m])
                            el = el1*el2
                            nominators_p.append([])
                            nominators_p[-1].append(expn*el)
                            nominators_p[-1].append(expm*el)
                            denominators_p.append(e[n]-e[m])
            dynamicObservable.lehmannNominators.update({statePair: allgather_list(nominators_p)})
            dynamicObservable.lehmannDenominators.update({statePair: allgather_list(denominators_p)})
        dynamicObservable.partitionFunction = self.getPartitionFunction()

    def setZeroFrequencyTerms(self, dynamicObservable):
        for statePair, operatorPair in dynamicObservable.operatorPairs.items():
            operator1 = operatorPair[0]
            operator2 = operatorPair[1]
            if type(operator1) != ndarray:
                operator1 = operator1.toarray()
                operator2 = operator2.toarray()
            e = self.getEnergyEigenvalues()
            ve = self.getEnergyEigenstates().toarray()
            fockstates = range(self.fockspaceSize)

            n_inds_p = array(scatter_list(fockstates))
            zeroFrequencyTerms_p = []
            for n in n_inds_p:
                n_op1 = ve[:,n].dot(operator1)
                op2_n = operator2.dot(ve[:,n])
                expn = exp(-self.beta*e[n])
                for m in fockstates:
                    if equals(e[n], e[m]):
                        el1 = n_op1.dot(ve[:,m])
                        if not equals(el1, 0):
                            el2 = ve[:,m].dot(op2_n)
                            if not equals(el2, 0):
                                el = el1*el2
                                zeroFrequencyTerms_p.append(-self.beta*expn*el)
            dynamicObservable.zeroFrequencyTerms.update({statePair: allgather_list(zeroFrequencyTerms_p)})

    def setLehmannTermsDynamicBosonic(self, dynamicObservable):
        self.setLehmannTermsDynamicFermionic(dynamicObservable)
        self.setZeroFrequencyTerms(dynamicObservable)

    def calcG1(self, singleParticleStatePairs = None):
        c = AnnihilationOperator(self.singleParticleBasis)
        if singleParticleStatePairs == None:
            statePairs = [(state1, state2) for state1, state2 in product(self.orderedSingleParticleStates, self.orderedSingleParticleStates)]
        else:
            statePairs = singleParticleStatePairs
        report('Calculating one-particle Green\'s function transition elements and energies...', self.verbose)
        t0 = time()
        for statePair in statePairs:
            c_state = c[statePair[0]]
            c_state_dag = c[statePair[1]].H
            self.g1.addOperatorPair(statePair, (c_state, c_state_dag))
        self.setLehmannTermsDynamic(self.g1)
        report('took '+str(time()-t0)[:4]+' seconds', self.verbose)

class GrandcanonicalEnsemble(CanonicalEnsemble):
    """use canonical with mu in H"""
    def __init__(self, hamiltonian, beta, mu, verbose = False):
        self.mu = mu
        CanonicalEnsemble.__init__(self, hamiltonian, beta, verbose)
        c = AnnihilationOperator(self.singleParticleBasis)
        muMatrix = mu * nsum([c[orb].H.dot(c[orb]) for orb in self.orderedSingleParticleStates], axis = 0)
        self.hamiltonian.matrix = self.hamiltonian.matrix - muMatrix
        self.filling = None

    def setMu(self, mu):
        c = AnnihilationOperator(self.singleParticleBasis)
        nMatrix = nsum([c[orb].H.dot(c[orb]) for orb in self.orderedSingleParticleStates], axis = 0)
        self.hamiltonian.matrix = self.hamiltonian.matrix + self.mu * nMatrix
        self.mu = mu
        self.hamiltonian.matrix = self.hamiltonian.matrix - mu * nMatrix
        self.energyEigenvalues = None
        self.energyEigenstates = None
        self.partitionFunction = None
        self.occupation = dict()
        report('Chemical potential set to '+str(mu), self.verbose)

    def setMuByFilling(self, filling, muMin, muMax, muTol = .001, maxiter = 100):
        c = AnnihilationOperator(self.singleParticleBasis)
        nMatrix = nsum([c[orb].H.dot(c[orb]) for orb in self.orderedSingleParticleStates], axis = 0)
        self.hamiltonian.matrix = self.hamiltonian.matrix + self.mu * nMatrix
        def fillingFunction(muTrial):
            self.hamiltonian.matrix = self.hamiltonian.matrix - muTrial * nMatrix
            self.calcEigensystem()
            self.calcPartitionFunction()
            self.calcOccupation()
            fillingTrial = self.getTotalOccupation()
            self.hamiltonian.matrix = self.hamiltonian.matrix + muTrial * nMatrix
            report('Filling(mu='+str(muTrial)+') = '+str(fillingTrial), self.verbose)
            self.filling = fillingTrial
            return fillingTrial - filling
        mu = bisect(fillingFunction, muMin, muMax, xtol = muTol, maxiter = maxiter)
        self.mu = mu
        self.hamiltonian.matrix = self.hamiltonian.matrix - mu * nMatrix
        report('Chemical potential set to '+str(mu), self.verbose)

class MicrocanonicalEnsemble(CanonicalEnsemble):
    def __init__(self, hamiltonian):
        CanonicalEnsemble.__init__(self, hamiltonian, 0)

    def calcPartitionFunction(self):
        self.partitionFunction = self.fockspaceSize
