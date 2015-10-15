from EasyED.ensembles import MicrocanonicalEnsemble, CanonicalEnsemble, GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard

import unittest, numpy

class TestEnsembles(unittest.TestCase):

    def runMicrocanonicalHubbardDimer(self):
        h = Hubbard([[0, -1], [-1, 0]], 1)
        ensemble = MicrocanonicalEnsemble(h)
        ensemble.calcOccupation()
        self.assertEqual(ensemble.getTotalOccupation(), 2)

    def runCanonicalHubbardDimer(self):
        h = Hubbard([[0, -1], [-1, 0]], 1)
        ensemble = CanonicalEnsemble(h, 1)
        ensemble.calcOccupation()
        self.assertTrue(ensemble.getTotalOccupation() < 2)
        
    def runGrandcanonicalHubbardDimer(self):
        h = Hubbard([[0, -1], [-1, 0]], 1)
        ensemble = GrandcanonicalEnsemble(h, 1, .5)
        ensemble.calcOccupation()
        self.assertEqual(ensemble.getTotalOccupation(), 2)
        for state in h.orderedSingleParticleStates:
            self.assertAlmostEqual(ensemble.occupation.getExpectationValue(state), .5)

    def runMuByFillingEstimation(self):
        h = Hubbard([[0, -1], [-1, 0]], 1)
        ensemble = GrandcanonicalEnsemble(h, 1, 0)
        ensemble.setMuByFilling(2, .4, .6)
        self.assertEqual(ensemble.getTotalOccupation(), 2)
        self.assertEqual(ensemble.mu, .5)

    def runG1Calculation(self):
        h = Hubbard([[0, -1], [-1, 0]], 1)
        ensemble = GrandcanonicalEnsemble(h, 1, .5)
        statePair = (('up', 1),('up', 1))
        ensemble.calcG1([statePair])
        ensemble.g1.setMesh(100, -2, 2)
        ensemble.g1.setRetarded(numpy.pi)
        a = ensemble.g1.getSpectralFunction(statePair)
        ensemble.g1.setMatsubaraMesh(20, 1)
        ensemble.g1.setMatsubara()
        giw = ensemble.g1.getMatsubara(statePair)
        #TODO need numerical data to benchmark
