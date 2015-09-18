from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load

us = [.75,1,3,8]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]

for u, fname in zip(us, fnames):
    betas = load(fname)[0,:]
    mus = load(fname)[1,:]
    t = -1
    r = 0
    for mu, beta in zip(mus, betas):
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.calcOccupation()
        report('N(beta='+str(beta)+',mu='+str(mu)+') = '+str(tetrahedron.getTotalOccupation()))
        report(str(tetrahedron.occupation.expectationValue.values()))
