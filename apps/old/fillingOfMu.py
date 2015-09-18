from numpy import linspace, save
from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard

beta = 100
t = -1
r = -1

us = [0, 1, 2]
datas = []
for u in us:
    h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = True)
    tetrahedron = GrandcanonicalEnsemble(h, beta, 0)
    datas.append([])
    mus = linspace(-2, 4, 60)
    fillings = list()
    datas[-1].append(mus)
    recalc = True
    for mu in mus:
        if recalc:
            tetrahedron.setMu(mu)
            tetrahedron.calcOccupation()
            n = tetrahedron.getTotalOccupation()
            if n < .001: recalc = False
        fillings.append(n)
    datas[-1].append(fillings)
fname = 'fillingOfMu.npy'
save(fname, datas)
print fname, 'ready'
