from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array

us = [1, 2.75, 3, 3.25, 4]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]

t = -1
r = -1
h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], 0, verbose = False)
c = AnnihilationOperator(h.singleParticleBasis)
sPlus0 = c['up', 0].H.dot(c['dn', 0])
sMinus0 = c['dn', 0].H.dot(c['up', 0])
sMinus1 = c['dn', 1].H.dot(c['up', 1])
chi_pm_hat = StaticObservable({'loc': sPlus0.dot(sMinus0), 'nn': sPlus0.dot(sMinus1)}, verbose = False)
results = list()

for u, fname in zip(us, fnames):
    betas = load(fname)[0,:]
    mus = load(fname)[1,:]
    chis_pm_loc = list()
    chis_pm_nn = list()
    for beta, mu in zip(betas, mus):
        report('u = '+str(u)+'; beta = '+str(beta)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannSumStatic(chi_pm_hat)
        chis_pm_loc.append(chi_pm_hat.getExpectationValue('loc'))
        chis_pm_nn.append(chi_pm_hat.getExpectationValue('nn'))
    results.append([betas, mus, chis_pm_loc, chis_pm_nn])

save('chi_temperature.npy', array(results))