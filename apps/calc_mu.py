from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import save, array, linspace

betas = linspace(1, 500, 500)
us = [0,0.75,1]
mu0 = 0
t = -1
r = -1

for u in us:
    h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
    mus = list()
    for beta in betas:
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu0, verbose = False)
        fnameMu = 'beta_mu_u'+str(u)+'.npy'
        tetrahedron.setMuByFilling(4, -10, 10)
        mus.append(tetrahedron.mu)
    save(fnameMu, array([betas, mus]))
