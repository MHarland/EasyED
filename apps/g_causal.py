from matplotlib import pyplot as plt
from tetrahedron import tetrahedron

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
tetrahedron.g1.setMesh(1000, -2, 2)

statePairs = [(('up', 0),('up', 0)), (('up', 0),('up', 1))]
for statePair in statePairs:
    tetrahedron.calcG1([statePair])
    tetrahedron.g1.setCausal()
    ax.plot(tetrahedron.g1.getMesh(), -tetrahedron.g1.getCausal(statePair).imag, label = str(statePair))

plt.legend()
plt.tight_layout()
plt.savefig('minus_im_g_causal.pdf', dpi=300)
plt.close()
