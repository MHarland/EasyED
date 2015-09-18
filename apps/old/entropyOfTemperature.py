from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.util import report
from numpy import load, save, array, log, linspace
from os.path import isfile

u = .75
t = -1
r = -1
h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
energies = list()
#temperatures = linspace(0.001, 1, 100)
#betas = list()
betas = linspace(1, 300, 300)
mus = list()
entropies = list()
zs = list()
#for i, temperature in enumerate(temperatures):
for i, beta in enumerate(betas):
    report(str(i)+'...')
    #beta = 1/temperature
    tetrahedron = GrandcanonicalEnsemble(h, beta, 0, verbose = False)
    tetrahedron.setMuByFilling(4, -1, 2)
    grand_h_hat = StaticObservable({'total': tetrahedron.hamiltonian.matrix}, False)
    tetrahedron.setLehmannSumStatic(grand_h_hat)
    avEnergy = grand_h_hat.getExpectationValue('total') - tetrahedron.hamiltonian.energyShift # artefact correction
    s = log(tetrahedron.getPartitionFunction()) + beta * avEnergy
    entropies.append(s)
    energies.append(avEnergy)
    zs.append(tetrahedron.getPartitionFunction())
    mus.append(tetrahedron.mu)
"""
for entropy,temperature,energy,z in zip(entropies, temperatures, energies, zs):
    report('', False)
    report('S(T='+str(temperature)+') = '+str(entropy), False)
    report('E(T='+str(temperature)+') = '+str(energy), False)
    report('Z(T='+str(temperature)+') = '+str(z), False)
"""
save('entropyOfTemperature.npy', array([betas, entropies, energies, zs, mus]))
save('beta_mu_u075.npy', array([betas, mus]))
