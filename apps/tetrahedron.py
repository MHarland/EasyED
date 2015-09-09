from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from EasyED.util import report
from numpy import load, save, array
from os.path import isfile

beta = 100
u = 0.75
mu0 = 0
t = -1
r = -1
h = Hubbard([[0,t,t,r],[t,0,r,t],[t,r,0,t],[r,t,t,0]], u, verbose = True)
tetrahedron = GrandcanonicalEnsemble(h, beta, mu0)

fnameMu = 'mu_u'+str(u)+'_b'+str(beta)+'.npy'
if isfile(fnameMu):
    report('loading chemical potential...')
    mu = load(fnameMu)[0]
    tetrahedron.setMu(mu)
    tetrahedron.calcEigensystem()
    tetrahedron.calcPartitionFunction()
else:
    tetrahedron.setMuByFilling(4, 0, 10)
    save(fnameMu, array([tetrahedron.mu]))
