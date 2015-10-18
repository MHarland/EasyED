from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import DynamicObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array, where, pi


us = [2.75, 3, 3.25]
us = [1,3,4]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
betas = array([10, 25, 44, 67, 100, 200])
betas = array([10, 100,200,400])
betas = array([5, 10, 15])
results = list()

for mu in mus:
    t = -1
    r = -1
    beta_mu = load(fname)
    results_u = list()
    for beta in betas:
        report('u = '+str(u)+'; beta = '+str(beta)+'; mu = '+str(mu)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        c = AnnihilationOperator(h.singleParticleBasis)
        sPlus0 = c['up', 0].H.dot(c['dn', 0])
        sMinus0 = c['dn', 0].H.dot(c['up', 0])
        sMinus1 = c['dn', 1].H.dot(c['up', 1])
        chi_pm = DynamicObservable({'loc': (sPlus0, sMinus0), 'nn': (sPlus0, sMinus1)}, False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannTermsDynamic(chi_pm)
        chi_pm.setMesh(1000, -4, 4)
        lehmannParams = [[1,-1],[1,1],pi/beta]
        results_mu.append([chi_pm.getMesh(), chi_pm.getCustom('loc', *lehmannParams).imag, chi_pm.getCustom('loc', *lehmannParams).real, chi_pm.getCustom('nn', *lehmannParams).imag, chi_pm.getCustom('nn', *lehmannParams).real])
    results.append(results_mu)
save('chi_pm.npy', array(results))
