from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import save, array, linspace

temperatures = linspace(.01, 3.01, 60)
betas = 1/temperatures
us = [1,2,3,4]
mu0 = 0
t = -1
r = -1

for u in us:
    mus = list()
    for beta in betas:
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu0, verbose = False)
        fnameMu = 'beta_mu_u'+str(u)+'b.npy'
        tetrahedron.setMuByFilling(4, -10, 30)
        mus.append(tetrahedron.mu)
        report('u = '+str(u)+', beta = '+str(beta)+' -> mu = '+str(tetrahedron.mu))
    save(fnameMu, array([betas, mus]))
