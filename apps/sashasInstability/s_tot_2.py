from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.observables import DynamicObservable
from EasyED.operators import AnnihilationOperator
from EasyED.util import report
from numpy import load, save, array, where, pi, sum as nsum

mus = [.27]
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
        stot = .5 * nsum([c['up', i].H.dot(c['up', i]) - c['dn', i].H.dot(c['dn', i]) for i in range(4)], axis = 0)
        s_2 = DynamicObservable({'tot': (stot, stot)}, False)
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.setLehmannTermsDynamic(s_2)
        s_2.setMesh(500, 0, 4)
        lehmannParams = [[1,-1],[1,1],pi/beta]
        results_mu.append([s_2.getMesh(), s_2.getCustom('tot', *lehmannParams).imag, s_2.getCustom('tot', *lehmannParams).real])
    results.append(results_mu)
save('s_tot_2.npy', array(results))
