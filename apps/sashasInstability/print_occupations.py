from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load

us = [3]
mus = [.26]
betas = [10]

for u in us:
    t = -1
    r = .3
    for mu, beta in zip(mus, betas):
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.calcOccupation()
        report('N(beta='+str(beta)+',mu='+str(mu)+') = '+str(tetrahedron.getTotalOccupation()))
        report(str(tetrahedron.occupation.expectationValue.values()))
