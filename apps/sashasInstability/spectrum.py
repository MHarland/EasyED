from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load, save, array, where, linspace

mus = [.27]
betas = array([20])
results = list()
u = 3

for mu in mus:
    t = -1
    r = .3
    results_u = list()
    for beta in betas:
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
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
