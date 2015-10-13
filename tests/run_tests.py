from tblocks import TestBlocks
from toperators import TestOperators
from thamiltonians import TestHamiltonians

import unittest

suite = unittest.TestSuite()
suite.addTest(TestBlocks('runBlockMatrix'))
suite.addTest(TestBlocks('runBlockMatrixProduct'))
suite.addTest(TestBlocks('runBlockState'))
suite.addTest(TestOperators('runSingleParticleBasis'))
suite.addTest(TestHamiltonians('runSpectrumComparisonWithSasha'))
suite.addTest(TestHamiltonians('runHubbardAtom'))
print
unittest.TextTestRunner(verbosity=2).run(suite)
