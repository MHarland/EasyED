from EasyED.hamiltonians import Hubbard
from EasyED.operators import AnnihilationOperator
from EasyED.util import report

import numpy
import unittest

class TestHamiltonians(unittest.TestCase):

    def runSpectrumComparisonWithSasha(self):
        almostZero = 10**(-10)
        mu = .27
        beta = 20
        results = list()
        u = 3
        t = -1
        r = .3
        h = Hubbard([[-mu,t,t,r],[t,-mu,r,t],[t,r,-mu,t],[r,t,t,-mu]], u, verbose = False)
        h.solve()
        energies = list(h.eigenEnergies)
        energiesSasha = numpy.loadtxt('spectrumSasha.dat')[:,1]
        self.assertEqual(len(energies), len(energiesSasha))
        for j, eS in enumerate(energiesSasha):
            energyFound = False
            for i, e in enumerate(energies):
                if abs(eS - e) < almostZero:
                    energyFound = True
                    del energies[i]
                    break
            self.assertTrue(energyFound)

    def runHubbardAtom(self):
        h = Hubbard([[-.5]], 1)
        h.solve()
        self.assertEqual(-.5, h.getGroundStateEnergy())
        self.assertEqual(h.getEnergies(), [0, .5])
        self.assertTrue((h.eigenEnergies == numpy.array([.5, 0, 0, .5])).all())
        for state in h.getGroundStatesAlgebraically():
            self.assertTrue(state in ['+1.0 c^(\'up\', 0)\n', '+1.0 c^(\'dn\', 0)\n'])
        for states_e, e in zip(h.getStatesEnergySortedAlgebraically(), h.getEnergies()):
            for state in states_e:
                if e == 0:
                    self.assertTrue(state in ['+1.0 c^(\'up\', 0)\n', '+1.0 c^(\'dn\', 0)\n'])
                elif e == .5:
                    self.assertTrue(state in ['+1.0 \n', '+1.0 c^(\'up\', 0) c^(\'dn\', 0)\n'])
                else:
                    self.assertTrue(False)
        
        c = AnnihilationOperator(h.getSingleParticleBasis())
        n_tot_hat = numpy.sum([c[s,0].H.dot(c[s,0]) for s in ['up', 'dn']], axis = 0)
        self.assertEqual(h.getGroundStates()[0].getQuantumNumber(n_tot_hat), 1)
        self.assertEqual(h.getGroundStates()[1].getQuantumNumber(n_tot_hat), 1)
