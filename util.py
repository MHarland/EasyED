from mpi4py import MPI as mpi
from numpy import array, dot as ndot, identity, matrix, sum as nsum

def diracDelta(x, y):
    return int(x == y)

def dot(*xs):
    xs = list(xs)
    res = array(xs.pop())
    for x in xs[::-1]:
        res = ndot(array(x), res)
    return res

def scatter_list(x):
    comm = mpi.COMM_WORLD
    cake = list()
    n_el = len(x)
    n_cpu = comm.size
    if comm.rank == 0:
        if n_el % n_cpu == 0:
            pieceSize = n_el/n_cpu
        else:
            pieceSize = n_el/n_cpu + 1
        for i, xi in enumerate(x):
            if i % pieceSize == 0:
                cake.append(list())
            cake[-1].append(xi)
    while len(cake) < n_cpu:
        cake.append(list())
    return comm.scatter(cake, root = 0)

def allgather_list(x):
    comm = mpi.COMM_WORLD
    cake = comm.allgather(x)
    y = list()
    for piece in cake:
        for xi in piece:
            y.append(xi)
    return y

def mpiSum(x):
    return mpi.COMM_WORLD.allreduce(x)

def sumScatteredLists(x):
    return nsum(allgather_list(x))
    #return mpiSum(x)

def report(s, verbose = True):
    if mpi.COMM_WORLD.Get_rank() == 0 and verbose:
        print s

def contains(numberList, number, accuracy):
    for l in numberList:
        if abs(number - l) < accuracy:
            return True
    return False

def getIndex(numberList, number, accuracy):
    for i, l in enumerate(numberList):
        if abs(number - l) < accuracy:
            return i
    return None
        
