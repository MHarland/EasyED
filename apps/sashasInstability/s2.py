from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import DynamicObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array, where, pi

u = 3
mus = [.26,.27,.28]
betas = array([20])
results = list()

for mu in mus:
    t = -1
    r = .3
    results_mu = list()
    for beta in betas:
        report('u = '+str(u)+'; beta = '+str(beta)+'; mu = '+str(mu)+'...')
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        c = AnnihilationOperator(h.singleParticleBasis)
        sz0 = (c['up', 0].H.dot(c['up', 0]) - c['dn', 0].H.dot(c['dn', 0])) * .5
        sz1 = (c['up', 1].H.dot(c['up', 1]) - c['dn', 1].H.dot(c['dn', 1])) * .5
        sz3 = (c['up', 3].H.dot(c['up', 3]) - c['dn', 3].H.dot(c['dn', 3])) * .5
        s2 = DynamicObservable({'loc': (sz0, sz0), 'nn': (sz0, sz1), 'nnn': (sz0, sz3)}, False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannTermsDynamic(s2)
        s2.setMesh(100, 0, 4)
        lehmannParams = [[1,-1],[1,1],pi/beta]
        results_mu.append([s2.getMesh(), s2.getCustom('loc', *lehmannParams).imag, s2.getCustom('loc', *lehmannParams).real, s2.getCustom('nn', *lehmannParams).imag, s2.getCustom('nn', *lehmannParams).real, s2.getCustom('nnn', *lehmannParams).imag, s2.getCustom('nnn', *lehmannParams).real])
    results.append(results_mu)
save('s2.npy', array(results))
