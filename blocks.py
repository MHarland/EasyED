from itertools import izip
from numpy import empty, dot as ndot

class BlockMatrix(object):
    def __init__(self, blocksizes):
        self.datablocks = list()
        self.blocksizes = blocksizes
        for size in blocksizes:
            self.datablocks.append(empty([size, size]))

    def __iter__(self):
        for block in self.datablocks:
            yield block

class BlockState(BlockMatrix):
    def __init__(self, blocksizes):
        self.datablocks = list()
        self.blocksizes = blocksizes
        for size in blocksizes:
            self.datablocks.append(empty([size]))    

def bdot(x, y):
    result = BlockMatrix(x.blocksizes)
    i = 0
    for xblock, yblock in izip(x, y):
        result.datablocks[i] = ndot(xblock, yblock)
        i += 1
    return result
