from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import DynamicObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array, where, pi

mus = [.26, .27, .28]
betas = array([20])
u = 3
t = -1
r = .3
results = list()

for mu in mus:
    results_mu = list()
    for beta in betas:
        report('u = '+str(u)+'; beta = '+str(beta)+'; mu = '+str(mu)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        c = AnnihilationOperator(h.singleParticleBasis)
        sz0 = .5 * (c['up', 0].H.dot(c['up', 0]) - c['dn', 0].H.dot(c['dn', 0]))
        sz1 = .5 * (c['up', 1].H.dot(c['up', 1]) - c['dn', 1].H.dot(c['dn', 1]))
        sz3 = .5 * (c['up', 3].H.dot(c['up', 3]) - c['dn', 3].H.dot(c['dn', 3]))
        chi_zz = DynamicObservable({'loc': (sz0, sz0), 'nn': (sz0, sz1), 'nnn': (sz0, sz3)}, 'bosonic')
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannTermsDynamic(chi_zz)
        chi_zz.setMesh(200, 0, 4)
        lehmannParams = [pi/beta,[1,-1]]
        results_mu.append([chi_zz.getMesh(), chi_zz.getCustom('loc', *lehmannParams).imag, chi_zz.getCustom('loc', *lehmannParams).real, chi_zz.getCustom('nn', *lehmannParams).imag, chi_zz.getCustom('nn', *lehmannParams).real, chi_zz.getCustom('nnn', *lehmannParams).imag, chi_zz.getCustom('nnn', *lehmannParams).real])
    results.append(results_mu)
save('chi_zz.npy', array(results))
