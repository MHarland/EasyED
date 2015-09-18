from matplotlib import pyplot as plt
from tetrahedron import tetrahedron

fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
tetrahedron.g1.setMatsubaraMesh(20, tetrahedron.beta)

statePairs = [(('up', 0),('up', 0)), (('up', 0),('up', 1))]
for statePair in statePairs:
    tetrahedron.calcG1([statePair])
    tetrahedron.g1.setMatsubara()
    ax1.plot(tetrahedron.g1.getMatsubaraMesh().imag, tetrahedron.g1.getMatsubara(statePair).imag, label = str(statePair))
    ax2.plot(tetrahedron.g1.getMatsubaraMesh().imag, tetrahedron.g1.getMatsubara(statePair).real, label = str(statePair))

ax1.legend()
ax1.set_ylabel('$\mathrm{Im}G(i\omega_n)$')
ax2.set_ylabel('$\mathrm{Re}G(i\omega_n)$')
ax2.set_xlabel('$\omega_n$')
plt.tight_layout()
plt.savefig('g_matsubara.pdf', dpi=300)
plt.close()
