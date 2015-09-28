from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import save, array, linspace

betas = linspace(1,400,400)
us = [3]
mu0 = 0
t = -1
r = .3

for u in us:
    mus = list()
    for beta in betas:
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu0, verbose = False)
        fnameMu = 'beta_mu_u'+str(u)+'.npy'
        tetrahedron.setMuByFilling(3, -10, 30)
        mus.append(tetrahedron.mu)
        report('u = '+str(u)+', beta = '+str(beta)+' -> mu = '+str(tetrahedron.mu))
    save(fnameMu, array([betas, mus]))
