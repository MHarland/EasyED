from itertools import izip, product
from numpy import zeros, dot as ndot

class Blocks(object):
    def __init__(self, blocksizes):
        self.datablocks = list()
        self.blocksizes = blocksizes
        for size in blocksizes:
            self.datablocks.append(zeros([size, size]))

    def __iter__(self):
        for block in self.datablocks:
            yield block

class BlockMatrix(Blocks):
    def setFullMatrixEntry(self, i, j, x):
        blockNr = 0
        blockrow = 0
        blockcol = 0
        origin = 0
        temp = 1
        while temp < i:
            blockNr += 1
            origin += self.blocksizes[blockNr]
            temp += self.blocksizes[blockNr+1]
        while origin + blockrow < i:
            blockrow += 1
        while origin + blockcol < j:
            blockcol += 1
        self.datablocks[blockNr][blockrow, blockcol] = x

class BlockState(Blocks):
    def __init__(self, blocksizes):
        self.datablocks = list()
        self.blocksizes = blocksizes
        for size in blocksizes:
            self.datablocks.append(zeros([size]))    

def bdot(x, y):
    result = BlockMatrix(x.blocksizes)
    i = 0
    for xblock, yblock in izip(x, y):
        result.datablocks[i] = ndot(xblock, yblock)
        i += 1
    return result
