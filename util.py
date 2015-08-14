from numpy import dot as ndot, identity, matrix, trace as ntrace

def diracDelta(x, y):
    return int(x == y)

def dot(*xs):
    res = identity(len(xs[0]))
    for x in xs:
        res = ndot(res, x)
    return matrix(res)

def trace(x):
    return ntrace(x) / float(len(x))

def scalar(*xs):
    return dot(*xs).flat[0]
