from tblocks import TestBlocks
from toperators import TestOperators

import unittest

suite = unittest.TestSuite()
suite.addTest(TestBlocks('runBlockMatrix'))
suite.addTest(TestBlocks('runBlockMatrixProduct'))
suite.addTest(TestBlocks('runBlockState'))
suite.addTest(TestOperators('runSingleParticleBasis'))
unittest.TextTestRunner(verbosity=2).run(suite)
