from EasyED.operators import AnnihilationOperator
from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load, save, array, where, linspace, sum as nsum

us = [3]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
betas = array([10])
results = list()
small = 10**(-10)

for u, fname in zip(us, fnames):
    t = -1
    r = -1
    beta_mu = load(fname)
    results_u = list()
    for beta in betas:
        h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = False)
        mu = beta_mu[1,where(beta_mu[0,:] == beta)[0][0]]
        report('u = '+str(u)+'; beta = '+str(beta)+'; mu = '+str(mu)+'...')
        tetrahedron = GrandcanonicalEnsemble(h, beta, mu, verbose = False)
        tetrahedron.calcEigensystem()
        c = AnnihilationOperator(tetrahedron.hamiltonian.getSingleParticleBasis())
        s_tot = list()
        for i in range(4):
            s_i = 3 * .5 * (c['up', i].H.dot(c['up', i]) - c['dn', i].H.dot(c['dn', i]))
            s_tot.append(s_i)
        s_tot = nsum(s_tot, axis = 0)
        s2_tot = s_tot.dot(s_tot)
        energies, degeneracies = tetrahedron.hamiltonian.getSpectrumEnergySorted()
        s2s = []
        datapointsDegeneracy = []
        datapoints = []
        for i, e, deg in zip(range(len(energies)), energies, degeneracies):
            n_states_e = len(tetrahedron.hamiltonian.getSuperpositionStatesEnergySorted()[i])
            s2s_i = [tetrahedron.hamiltonian.getSuperpositionStatesEnergySorted()[i][j].getQuantumNumber(s2_tot) for j in range(n_states_e)]
            for s2 in s2s_i:
                is_element = False
                for datapoint in datapoints:
                    if abs(s2-datapoint[0]) < small and abs(e-datapoint[1]) < small:
                        is_element = True
                        break
                if not is_element:
                    s2s.append(s2)
                    datapoints.append((s2, e))
                    datapointsDegeneracy.append(1)
                else:
                    k = 0
                    for k, datapoint in enumerate(datapoints):
                        if abs(s2-datapoint[0]) < small and abs(e-datapoint[1]) < small:
                            break
                    datapointsDegeneracy[k] += 1
        datapoints = array(datapoints)
        results_u.append([datapoints[:,0], datapoints[:,1], datapointsDegeneracy])
    results.append(results_u)
save('energy_s2.npy', array(results))
