from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load, save, array, where, linspace

us = [2.75,3,3.25]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
betas = array([10, 15, 20, 25, 30])
#betas = array([10,69,70,86,87])
results = list()

for u, fname in zip(us, fnames):
    t = -1
    r = -1
    beta_mu = load(fname)
    results_u = list()
    for beta in betas:
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        mu = beta_mu[1,where(beta_mu[0,:] == beta)[0][0]]
        report('u = '+str(u)+'; beta = '+str(beta)+'; mu = '+str(mu)+'...')
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.calcEigensystem()
        energies = list()
        degeneracies = list()
        for e, d in zip(*tetrahedron.hamiltonian.getSpectrum(10**(-10))):
            energies.append(e)
            degeneracies.append(d)
        results_u.append([energies, degeneracies])
    results.append(results_u)
save('spectrum.npy', array(results))
