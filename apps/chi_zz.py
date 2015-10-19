from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import DynamicObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array, where, pi

betas = array([10, 15, 20, 25, 30])
#betas = array([10,69,70,86,87])
us = [2.75,3,3.25]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
t = -1
r = -1
results = list()

for u, fname in zip(us, fnames):
    mu_beta = load(fname)
    results_u = list()
    for beta in betas:
        mu = mu_beta[1,where(mu_beta[0,:] == beta)[0][0]]
        report('u = '+str(u)+'; beta = '+str(beta)+'; mu = '+str(mu)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        c = AnnihilationOperator(h.singleParticleBasis)
        sz0 = .5 * (c['up', 0].H.dot(c['up', 0]) - c['dn', 0].H.dot(c['dn', 0]))
        sz1 = .5 * (c['up', 1].H.dot(c['up', 1]) - c['dn', 1].H.dot(c['dn', 1]))
        chi_zz = DynamicObservable({'loc': (sz0, sz0), 'nn': (sz0, sz1)}, 'bosonic')
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannTermsDynamic(chi_zz)
        chi_zz.setMesh(400, 0, 8)
        lehmannParams = [pi/beta,[1,-1]]
        results_u.append([chi_zz.getMesh(), chi_zz.getCustom('loc', *lehmannParams).imag, chi_zz.getCustom('loc', *lehmannParams).real, chi_zz.getCustom('nn', *lehmannParams).imag, chi_zz.getCustom('nn', *lehmannParams).real])
    results.append(results_u)
save('chi_zz.npy', array(results))
