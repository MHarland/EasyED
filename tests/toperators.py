from EasyED.operators import SingleParticleBasis

import unittest

class TestOperators(unittest.TestCase):
    """
    def __init__(self, *args, **kwargs):
        super(TestOperators, self).__init__(*args, **kwargs)
    """

    def runSingleParticleBasis(self):
        basis = SingleParticleBasis([['u', 'd'], [0, 1]])
        self.assertTrue(basis.getStateAlgebraically(13) == 'c^(\'u\', 0) c^(\'d\', 0) c^(\'d\', 1)')
        self.assertTrue(basis.getFockspaceNr((1,0,1,1)) == 13)
        self.assertTrue(basis.getOccupationRep(13) == '1011')
        self.assertTrue(basis.getFockspaceNr((0,1,0,0)) == 2)
        self.assertTrue(basis.getSingleParticleStateNr('u', 1) == 1)
        self.assertTrue(basis.getOccupationRep(2) == '0100')
        for sps, fockNr in zip(basis.orderedSingleParticleStates, [1,2,4,8]):
            self.assertTrue(basis.getFockspaceNr(singleParticleState = sps) == fockNr)

    def runAnnihilationOperator(self):
        pass
        
        
        
