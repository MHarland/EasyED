from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import StaticObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array

us = [1,2,3,4]
fnames = ['beta_mu_u'+str(u)+'b.npy' for u in us]

t = -1
r = -1
h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], 0, verbose = False)
c = AnnihilationOperator(h.singleParticleBasis)
sZ0 = .5 * (c['up', 0].H.dot(c['up', 0]) - c['dn', 0].H.dot(c['dn', 0]))
sZ1 = .5 * (c['up', 1].H.dot(c['up', 1]) - c['dn', 1].H.dot(c['dn', 1]))
chi_zz_hat = StaticObservable({'loc': sZ0.dot(sZ0), 'nn': sZ0.dot(sZ1)}, verbose = False)
results = list()

for u, fname in zip(us, fnames):
    betas = load(fname)[0,:]
    mus = load(fname)[1,:]
    s2_loc = list()
    s2_nn = list()
    for beta, mu in zip(betas, mus):
        report('u = '+str(u)+'; beta = '+str(beta)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannSumStatic(chi_zz_hat)
        s2_loc.append(3*chi_zz_hat.getExpectationValue('loc'))
        s2_nn.append(3*chi_zz_hat.getExpectationValue('nn'))
    results.append([betas, mus, s2_loc, s2_nn])

save('s2_temperature.npy', array(results))
