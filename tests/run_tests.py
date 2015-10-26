from tblocks import TestBlocks
from toperators import TestOperators
from thamiltonians import TestHamiltonians
from tensembles import TestEnsembles

import unittest

suite = unittest.TestSuite()
suite.addTest(TestBlocks('runBlockMatrix'))
suite.addTest(TestBlocks('runBlockMatrixProduct'))
suite.addTest(TestBlocks('runBlockState'))
suite.addTest(TestOperators('runSingleParticleBasis'))
suite.addTest(TestHamiltonians('runSpectrumComparisonWithSasha'))
suite.addTest(TestHamiltonians('runHubbardAtom'))
suite.addTest(TestHamiltonians('runEigenstatesEquation'))
suite.addTest(TestEnsembles('runMicrocanonicalHubbardDimer'))
suite.addTest(TestEnsembles('runCanonicalHubbardDimer'))
suite.addTest(TestEnsembles('runGrandcanonicalHubbardDimer'))
suite.addTest(TestEnsembles('runMuByFillingEstimation'))
suite.addTest(TestEnsembles('runG1Calculation'))
print
unittest.TextTestRunner(verbosity=2).run(suite)
