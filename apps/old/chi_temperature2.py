from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array, linspace
from tetrahedron import tetrahedron

u = .75
t = -1
r = -1
h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
temperatures = linspace(.001, .3, 300)
chis_pm_loc = list()
chis_pm_nn = list()
mus = list()

c = AnnihilationOperator(tetrahedron.singleParticleBasis)

sPlus0 = c['up', 0].H.dot(c['dn', 0])
sMinus0 = c['dn', 0].H.dot(c['up', 0])
sMinus1 = c['dn', 1].H.dot(c['up', 1])
chi_pm_hat = StaticObservable({'loc': sPlus0.dot(sMinus0), 'nn': sPlus0.dot(sMinus1)}, verbose = False)

for i, t in enumerate(temperatures):
    report(str(i)+'...')
    beta = 1./t
    tetrahedron = GrandcanonicalEnsemble(h, beta, 0, verbose = False)
    tetrahedron.setMuByFilling(4, -1, 2)
    tetrahedron.setLehmannSumStatic(chi_pm_hat)
    chis_pm_loc.append(chi_pm_hat.getExpectationValue('loc'))
    chis_pm_nn.append(chi_pm_hat.getExpectationValue('nn'))
    mus.append(tetrahedron.mu)

save('chi_temperature2.npy', array([temperatures, chis_pm_loc, chis_pm_nn]))
save('mu_temperature_u075b.npy', array([temperatures, mus]))
