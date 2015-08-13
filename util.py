from numpy import dot as ndot, identity

def diracDelta(x, y):
    return int(x == y)

def dot(*xs):
    res = identity(len(xs[0]))
    for x in xs:
        res = ndot(res, x)
    return res
