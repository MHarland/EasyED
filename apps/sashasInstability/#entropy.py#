from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.util import report
from numpy import load, save, array, log, linspace

mus = [.26,.27,.28]
temperatures = linspace(.0001,.1,20)
betas = 1/temperatures
u = 3
t = -1
r = .3
results = list()

for mu in mus:
    energies = list()
    entropies = list()
    zs = list()
    for beta in betas:
        report('mu = '+str(mu)+'; beta = '+str(beta)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        grand_h_hat = StaticObservable({'total': tetrahedron.hamiltonian.matrix}, False)
        tetrahedron.setLehmannSumStatic(grand_h_hat)
        avEnergy = grand_h_hat.getExpectationValue('total') - tetrahedron.hamiltonian.energyShift # artefact correction
        entropies.append(log(tetrahedron.getPartitionFunction()) + beta * avEnergy)
        energies.append(avEnergy)
        zs.append(tetrahedron.getPartitionFunction())
    results.append([betas, entropies, energies, zs])

save('entropy.npy', array(results))
