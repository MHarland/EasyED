from EasyED.hamiltonians import Hubbard
from EasyED.operators import AnnihilationOperator
from EasyED.util import report

from itertools import product
from matplotlib import pyplot as plt
from numpy import array, sqrt, sort, arange
import numpy



numpy.set_printoptions(suppress=True)
structure = Hubbard([[.5, 0], [0, .5]], 1)
c = AnnihilationOperator(structure.getSingleParticleBasis())

report(structure.blocksizes)
structure.solve()
print structure.eigenEnergies

print 'groundstateenergy: ', structure.getGroundStateEnergy()
print 'groundstates:'
for state in structure.getGroundStatesAlgebraically():
    print state
    print

print 'states:'
for energyGroup, energy in zip(structure.getStatesEnergySortedAlgebraically(), structure.getEnergies()):
    print 'E = ', energy
    print
    for state in energyGroup:
        print state
        print
    print

sz_tot = .5 * (c['up', 0].H.dot(c['up', 0]) - c['dn', 0].H.dot(c['dn', 0])) + .5 * (c['up', 1].H.dot(c['up', 1]) - c['dn', 1].H.dot(c['dn', 1]))
chi_zz_tot = sz_tot.dot(sz_tot)
print 'corresponding S2_tot qnrs:'
energies, degeneracies = structure.getSpectrum()
for energyGroup, energy  in zip(structure.getSuperpositionStatesEnergySorted(), sort(energies)):
    for state in energyGroup:
        print str(energy)+' -> '+str(3*state.getQuantumNumber(chi_zz_tot))
print
n_tot = c['up', 0].H.dot(c['up', 0]) + c['dn', 0].H.dot(c['dn', 0]) + c['up', 1].H.dot(c['up', 1]) + c['dn', 1].H.dot(c['dn', 1])
print 'corresponding n_tot qnrs:'
energies, degeneracies = structure.getSpectrum()
for energyGroup, energy  in zip(structure.getSuperpositionStatesEnergySorted(), sort(energies)):
    for state in energyGroup:
        print str(energy)+' -> '+str(state.getQuantumNumber(n_tot))

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
es = list()
for e, d in zip(*structure.getSpectrum()):
    es.append(e)
    ax.plot([e]*2, [0, d], color = 'blue')
ax.set_xlim(min(es)-1,max(es)+1)
plt.tight_layout()
plt.savefig('testSpectrum.pdf', dpi=300)
plt.close()

#sz_1 = .5 * (c['up', 0].H.dot(c['up', 0]) - c['dn', 0].H.dot(c['dn', 0]))
#sz_2 = .5 * (c['up', 1].H.dot(c['up', 1]) - c['dn', 1].H.dot(c['dn', 1]))
#chi_zz_tot = sz_1.dot(sz_1) + sz_2.dot(sz_2)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
energies, degeneracies = structure.getSpectrumEnergySorted()
s2s = []
datapointsDegeneracy = []
datapoints = []
for i, e, deg in zip(range(len(energies)), energies, degeneracies):
    for s2 in [3*structure.getSuperpositionStatesEnergySorted()[i][j].getQuantumNumber(chi_zz_tot) for j in range(len(structure.getSuperpositionStatesEnergySorted()[i]))]:
        if not (s2, e) in datapoints:
            ax.scatter(s2, e, marker = '_', color = 'black')
            s2s.append(s2)
            datapoints.append((s2, e))
            datapointsDegeneracy.append(1)
        else:
            k = 0
            for k, datapoint in enumerate(datapoints):
                if (s2, e) == datapoint:
                    break
            datapointsDegeneracy[k] += 1
for datapoint, deg in zip(datapoints, datapointsDegeneracy):
    ax.text(datapoint[0], datapoint[1], str(deg), color = 'gray')
ax.set_xticks(s2s)
ax.set_xlim(0-.1*max(s2s),1.1*max(s2s))
ax.set_ylim(0-.1*max(energies),1.1*max(energies))
ax.set_xlabel('$S^2_{tot}$')
ax.set_ylabel('$E$')
plt.tight_layout()
plt.savefig('test_edegs2.pdf', dpi=300)
plt.close()
