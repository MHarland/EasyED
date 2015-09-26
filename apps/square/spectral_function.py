from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load, save, array, where, pi

us = [1,8,12]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
betas = array([1,5,10,15])
results = list()

for u, fname in zip(us, fnames):
    t = -1
    r = 0
    beta_mu = load(fname)
    results_u = list()
    for beta in betas:
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        mu = beta_mu[1,where(beta_mu[0,:] == beta)[0][0]]
        report('u = '+str(u)+'; beta = '+str(beta)+'; mu = '+str(mu)+'...')
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.g1.setMesh(200,-10,10)
        tetrahedron.calcG1([(('up',0),('up',0))])
        tetrahedron.g1.setRetarded(pi/beta)
        results_u.append([tetrahedron.g1.getMesh(), tetrahedron.g1.getSpectralFunction((('up',0),('up',0)))])
    results.append(results_u)
save('spectral_function.npy', array(results))
