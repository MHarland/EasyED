from itertools import product
from numpy import trace, exp, sum as nsum, identity, where
from operators import AnnihilationOperator
from util import dot, scatter_list, sumScatteredLists, report

class CanonicalEnsemble(object):
    def __init__(self, hamiltonian, beta):
        self.hamiltonian = hamiltonian
        self.singleParticleBasis = hamiltonian.singleParticleBasis
        self.orderedSingleParticleStates = hamiltonian.orderedSingleParticleStates
        self.fockspaceSize = hamiltonian.fockspaceSize
        self.partitionFunction = None
        self.beta = beta
        self.energyEigenstates = None
        self.energyEigenvalues = None
        self.occupation = dict()
        self.verbose = True

    def calcOccupation(self):
        c = AnnihilationOperator(self.singleParticleBasis)
        e = self.getEnergyEigenvalues()
        ve = self.getEnergyEigenstates()
        z = self.getPartitionFunction()
        inds = self.getIndsForSum()

        report('Calculating Occupation...', self.verbose)
        for state in self.orderedSingleParticleStates:
            inds_p = scatter_list(inds)
            terms_p = list()
            for m, n in inds_p:
                el = dot(ve[:,m],c[state].H,ve[:,n])
                terms_p.append(exp(-self.beta*e[n])*el*el.conjugate())
            self.occupation.update({state: sumScatteredLists(terms_p)/z})

    def getTotalOccupation(self):
        return nsum(self.occupation.values())

    def calcPartitionFunction(self):
        report('Calculating Partition Function...', self.verbose)
        self.partitionFunction = nsum([exp(-self.beta*eva) for eva in self.getEnergyEigenvalues()])

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
        report('Calculating eigensystem...', self.verbose)
        self.hamiltonian.solve()
        self.energyEigenstates = self.hamiltonian.eigenStates
        self.energyEigenvalues = self.hamiltonian.eigenEnergies

"""    use canonical with mu in H
class GrandcanonicalEnsemble(CanonicalEnsemble):
    def __init__(self, hamiltonian, beta, mu):
        CanonicalEnsemble.__init__(self, hamiltonian, beta)
        c = AnnihilationOperator(self.singleParticleBasis)
        self.hamiltonian.matrix -= mu * nsum([dot(c[orb].H, c[orb]) for orb in self.orderedSingleParticleStates], axis = 0)
"""

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
