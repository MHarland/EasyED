from itertools import izip, product
from numpy import zeros, dot as ndot, asmatrix, sum as nsum

class Blocks(object):
    def __init__(self, blocksizes):
        self.datablocks = list()
        self.blocksizes = blocksizes

    def __iter__(self):
        for block in self.datablocks:
            yield block

    def __getitem__(self, x):
        return self.datablocks[x]

    def show(self):
        for block in self.datablocks:
            print block

    def getBlockCoords(self, i, j = None):
        blockNr = 0
        blockrow = 0
        blockcol = 0
        origin = 0
        temp = self.blocksizes[0]
        while temp <= i:
            origin += self.blocksizes[blockNr]
            temp += self.blocksizes[blockNr+1]
            blockNr += 1
        while origin + blockrow < i:
            blockrow += 1
        if j != None:
            assert origin <= i < temp and origin <= j < temp, 'Cannot access values outside the blocks.'
            while origin + blockcol < j:
                blockcol += 1
            return blockNr, blockrow, blockcol
        else:
            return blockNr, blockrow

    def getBlock(self, blockNr):
        return self.datablocks[blockNr]

class BlockMatrix(Blocks):
    def __init__(self, blocksizes):
        Blocks.__init__(self, blocksizes)
        for size in blocksizes:
            self.datablocks.append(asmatrix(zeros([size, size])))
        self.shape = [nsum(blocksizes)]*2

    def __setitem__(self, key, item):
        blockNr, blockrow, blockcol = self.getBlockCoords(*key)
        self.datablocks[blockNr][blockrow, blockcol] = item

    def __getitem__(self, key):
        blockNr, blockrow, blockcol = self.getBlockCoords(*key)
        return self.datablocks[blockNr][blockrow, blockcol]

    def bdot(self, y):
        return bdot(self, y)

class BlockState(Blocks):
    def __init__(self, blocksizes):
        Blocks.__init__(self, blocksizes)
        for size in blocksizes:
            self.datablocks.append(zeros([size]))
        self.shape = [nsum(blocksizes)]

    def __setitem__(self, key, item):
        blockNr, blockrow = self.getBlockCoords(key)
        self.datablocks[blockNr][blockrow] = item

    def __getitem__(self, key):
        blockNr, blockrow = self.getBlockCoords(key)
        return self.datablocks[blockNr][blockrow]

    def bdot(self, y):
        return bdot(self, y)

def bdot(x, y):
    assert len(x.shape) < 3 and len(y.shape) < 3, 'wrong shape'
    if len(x.shape) == 2 and len(y.shape) == 2:
        result = BlockMatrix(x.blocksizes)
    elif len(x.shape) == len(y.shape) == 1:
        result = 0
        for i in range(x.shape[0]):
            result += x[i] * y[i]
        return result
    else:
        result = BlockState(x.blocksizes)
    i = 0
    for xblock, yblock in izip(x, y):
        result.datablocks[i] = ndot(xblock, yblock)
        i+= 1
    return result
