from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import DynamicObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array, where, pi

us = [1,8,12]
betas = array([1,5,10])
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
results = list()

for u, fname in zip(us, fnames):
    t = -1
    r = 0
    beta_mu = load(fname)
    results_u = list()
    for beta in betas:
        eta = pi/beta
        mu = beta_mu[1,where(beta_mu[0,:] == beta)[0][0]]
        report('u = '+str(u)+'; beta = '+str(beta)+'; mu = '+str(mu)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        c = AnnihilationOperator(h.singleParticleBasis)
        sPlus0 = c['up', 0].H.dot(c['dn', 0])
        sMinus0 = c['dn', 0].H.dot(c['up', 0])
        sMinus1 = c['dn', 1].H.dot(c['up', 1])
        chi_pm = DynamicObservable({'loc': (sPlus0, sMinus0), 'nn': (sPlus0, sMinus1)}, False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannTermsDynamic(chi_pm)
        chi_pm.setMesh(500, -10, 10)
        results_u.append([chi_pm.getMesh(), chi_pm.getCustom('loc', [1, -1], [1, 0], eta).imag, chi_pm.getCustom('loc', [1, -1], [1, 0], eta).real, chi_pm.getCustom('nn', [1, -1], [1, 0], eta).imag, chi_pm.getCustom('nn', [1, -1], [1, 0], eta).real])
    results.append(results_u)
save('chi_pm.npy', array(results))
