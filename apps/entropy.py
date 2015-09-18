from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.util import report
from numpy import load, save, array, log

us = [1,2.25,2.5,2.75,3,3.25,4]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
t = -1
r = -1
results = list()

for u, fname in zip(us, fnames):
    betas = load(fname)[0,:]
    mus = load(fname)[1,:]
    energies = list()
    entropies = list()
    zs = list()
    for beta, mu in zip(betas, mus):
        report('u = '+str(u)+'; beta = '+str(beta)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        grand_h_hat = StaticObservable({'total': tetrahedron.hamiltonian.matrix}, False)
        tetrahedron.setLehmannSumStatic(grand_h_hat)
        avEnergy = grand_h_hat.getExpectationValue('total') - tetrahedron.hamiltonian.energyShift # artefact correction
        entropies.append(log(tetrahedron.getPartitionFunction()) + beta * avEnergy)
        energies.append(avEnergy)
        zs.append(tetrahedron.getPartitionFunction())
    results.append([betas, mus, entropies, energies, zs])

save('entropy.npy', array(results))
