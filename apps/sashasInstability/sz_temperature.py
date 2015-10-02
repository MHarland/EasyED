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
c = AnnihilationOperator(h.singleParticleBasis)
sZ = .5 * (c['up', 0].H.dot(c['up', 0]) - c['dn', 0].H.dot(c['dn', 0]))
s_z_hat = StaticObservable({'loc': sZ}, verbose = False)
results = list()

for u, fname in zip(us, fnames):
    betas = load(fname)[0,:]
    mus = load(fname)[1,:]
    sz_loc = list()
    for beta, mu in zip(betas, mus):
        report('u = '+str(u)+'; beta = '+str(beta)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannSumStatic(s_z_hat)
        sz_loc.append(s_z_hat.getExpectationValue('loc'))
    results.append([betas, mus, sz_loc])

save('sz_temperature.npy', array(results))
