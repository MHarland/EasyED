from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import save, array, linspace

betas = array([1,5,10,15])
us = [1,4,8,12]
mu0 = 0
t = -1
r = 0

for u in us:
    mus = list()
    for beta in betas:
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu0, verbose = False)
        fnameMu = 'beta_mu_u'+str(u)+'.npy'
        tetrahedron.setMuByFilling(4, 0, u)
        mus.append(tetrahedron.mu)
        report('u = '+str(u)+', beta = '+str(beta)+' -> mu = '+str(tetrahedron.mu))
    save(fnameMu, array([betas, mus]))
