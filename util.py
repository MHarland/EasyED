from numpy import array, dot as ndot, identity, matrix, trace as ntrace

def diracDelta(x, y):
    return int(x == y)

def dot(*xs):
    xs = list(xs)
    res = array(xs.pop())
    for x in xs[::-1]:
        res = ndot(array(x), res)
    return res

def trace(x):
    return ntrace(x) / float(len(x))
