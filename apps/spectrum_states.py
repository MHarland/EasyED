from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load, save, array, where, linspace

us = [3]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
betas = array([10, 100, 200])
energyIndices_beta = [[],[],[0,1,2]]
results = list()

for u, fname in zip(us, fnames):
    t = -1
    r = -1
    beta_mu = load(fname)
    results_u = list()
    for beta, energyIndices in zip(betas, energyIndices_beta):
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        mu = beta_mu[1,where(beta_mu[0,:] == beta)[0][0]]
        report('u = '+str(u)+'; beta = '+str(beta)+'; mu = '+str(mu)+':')
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.calcEigensystem()
        for energyIndex in energyIndices:
            report('state nr. '+str(energyIndex)+' with E = ' +str(tetrahedron.hamiltonian.getEnergies()[energyIndex])+':\n')
            for state in tetrahedron.hamiltonian.getStatesEnergySortedAlgebraically()[energyIndex]:
                report(state)
                report('')
        report('')
