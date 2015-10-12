from EasyED.blocks import BlockMatrix, BlockState, bdot

from numpy import array
from itertools import product
import unittest

class TestBlocks(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestBlocks, self).__init__(*args, **kwargs)
        self.x_bm = BlockMatrix([2,2,1])
        a = array([[1,2],[3,4]])
        b = array([[1,5],[5,4]])
        c = array([[10]])
        abc = [a, b, c]
        for i, block in enumerate(self.x_bm):
            self.x_bm.datablocks[i] = abc[i]
        self.x_numpy = array([[1,2,0,0,0],[3,4,0,0,0],[0,0,1,5,0],[0,0,5,4,0],[0,0,0,0,10]])

    def runBlockMatrix(self):
        k = 0
        for i, j in product(range(5), range(5)):
            if self.x_numpy[i, j] != 0:
                k += 1
                self.assertTrue(self.x_numpy[i, j] == self.x_bm[i, j])
        self.assertTrue(k == 9)
        self.assertTrue((self.x_bm.getBlock(1) == array([[1,5],[5,4]])).all())
        self.assertTrue(self.x_bm[2,3] == 5)

    def runBlockMatrixProduct(self):
        xx_numpy = self.x_numpy.dot(self.x_numpy)
        xx_bm = bdot(self.x_bm, self.x_bm)
        for i, j in product(range(5), range(5)):
            if xx_numpy[i, j] != 0:
                self.assertTrue(xx_numpy[i, j] == xx_bm[i, j])

    def runBlockState(self):
        a = [1,2]
        b = [3,4]
        c = [5]
        d_numpy = array(a+b+c)
        d_bs = BlockState([2,2,1])
        for i, dataBlock in zip(range(3), [a,b,c]):
            d_bs.datablocks[i] = array(dataBlock)
        for i in range(5):
            self.assertTrue(d_numpy[i] == d_bs[i])
        d2_bs = BlockState([2,2,1])
        for i in range(5):
            d2_bs[i] = d_bs[i]
        for i in range(5):
            self.assertTrue(d_numpy[i] == d2_bs[i])
        xd_numpy = self.x_numpy.dot(d_numpy)
        xd_bms = bdot(self.x_bm, d_bs)
        for i in range(5):
            self.assertTrue(xd_numpy[i] == xd_bms[i])
        dd_numpy = d_numpy.dot(d_numpy)
        dd_bs = bdot(d_bs, d_bs)
        self.assertTrue(dd_numpy == dd_bs)

        
