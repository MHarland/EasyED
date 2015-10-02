from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array, linspace

mus = [.26,.27,.28]
temperatures = linspace(.0001,.1,20)
betas = 1/temperatures
u = 3
t = -1
r = .3
h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], 0, verbose = False)
results = list()

for mu in mus:
    specificHeat = list()
    for beta  in betas:
        report('mu = '+str(mu)+'; beta = '+str(beta)+'...')
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
