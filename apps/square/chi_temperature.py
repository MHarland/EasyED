from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array, linspace

us = [0.75, 1, 2, 4]
betas = linspace(1,400,400)

t = -1
r = 0
h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], 0, verbose = False)
c = AnnihilationOperator(h.singleParticleBasis)
sPlus0 = c['up', 0].H.dot(c['dn', 0])
sMinus0 = c['dn', 0].H.dot(c['up', 0])
sMinus1 = c['dn', 1].H.dot(c['up', 1])
chi_pm_hat = StaticObservable({'loc': sPlus0.dot(sMinus0), 'nn': sPlus0.dot(sMinus1)}, verbose = False)
results = list()

for u in us:
    mu = u * .5
    chis_pm_loc = list()
    chis_pm_nn = list()
    for beta in betas:
        report('u = '+str(u)+'; beta = '+str(beta)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannSumStatic(chi_pm_hat)
        chis_pm_loc.append(chi_pm_hat.getExpectationValue('loc'))
        chis_pm_nn.append(chi_pm_hat.getExpectationValue('nn'))
    results.append([betas, chis_pm_loc, chis_pm_nn])

save('chi_temperature.npy', array(results))
