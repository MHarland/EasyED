from itertools import product
from numpy import trace, exp, sum as nsum, identity
from operators import AnnihilationOperator
from util import scalar

class CanonicalEnsemble(object):
    def __init__(self, hamiltonian, beta):
        self.hamiltonian = hamiltonian
        self.singleParticleBasis = hamiltonian.singleParticleBasis
        self.partitionFunction = None
        self.beta = beta
        self.energyEigenstates = None
        self.energyEigenvalues = None

    def occupation(self, *orbital):
        c = AnnihilationOperator(self.singleParticleBasis)
        z = self.getPartitionFunction()
        e = self.getEnergyEigenvalues()
        ve = self.getEnergyEigenstates()
        inds = range(len(ve))
        return nsum([exp(-self.beta*e[n])*scalar(ve[n],c[orbital].H,ve[m])*scalar(ve[m],c[orbital],ve[n]) for m, n in product(inds, inds)]) / z

    def calcPartitionFunction(self):
        self.partitionFunction = nsum([exp(-self.beta*eva) for eva in self.getEnergyEigenvalues()])

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
    def __init__(self, hamiltonian, beta, mu):
        CanonicalEnsemble.__init__(self, hamiltonian, beta)
        self.hamiltonian.matrix -= mu * identity(len(self.hamiltonian.matrix))
