from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array

us = [3]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]

t = -1
r = .3
h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], 0, verbose = False)
results = list()

for u, fname in zip(us, fnames):
    betas = load(fname)[0,:]
    mus = load(fname)[1,:]
    specificHeat = list()
    for beta, mu in zip(betas, mus):
        report('u = '+str(u)+'; beta = '+str(beta)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        eMatrix = h.matrix.copy()
        eSquaredMatrix = eMatrix.dot(eMatrix)
        e_hat = StaticObservable({'total': eMatrix}, verbose = False)
        eSquared_hat = StaticObservable({'total': eSquaredMatrix}, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannSumStatic(e_hat)
        tetrahedron.setLehmannSumStatic(eSquared_hat)
        e = e_hat.getExpectationValue('total')
        eSquared = eSquared_hat.getExpectationValue('total')
        specificHeat.append(beta**2 * (eSquared - e**2))
    results.append([betas, specificHeat])

save('specific_heat_temperature.npy', array(results))
