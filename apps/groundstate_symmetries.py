from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load, save, array, where

us = [1, 2, 3, 4]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
betas = array([10])
results = list()

for u, fname in zip(us, fnames):
    t = -1
    r = -1
    beta_mu = load(fname)
    for beta in betas:
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        mu = beta_mu[1,where(beta_mu[0,:] == beta)[0][0]]
        report('u = '+str(u)+', beta = '+str(beta)+', mu = '+str(mu)+':')
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.calcEigensystem()
        for s in tetrahedron.hamiltonian.getGroundStates(threshold = .3):
            report(s)
            report('')
        report('')
        report('')