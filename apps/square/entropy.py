from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.util import report
from numpy import load, save, array, log, linspace

us = [.75, 1]
betas = linspace(1, 400, 400)
t = -1
r = 0
results = list()

for u in us:
    mu = u * .5
    energies = list()
    entropies = list()
    zs = list()
    for beta in betas:
        report('u = '+str(u)+'; beta = '+str(beta)+'...')
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
