from numpy import load, array
from matplotlib import pyplot as plt
from matplotlib_to_latex import set_poster_parameters as set_mpl

a_u_beta = load('spectral_function.npy')
us = [2.75, 3, 3.25]
us = [.75, 1, 2, 2.75, 3, 3.25, 4]
us = [1,3,4]
u_inds = range(3) # indices of us that will be plotted
betas = array([10, 100, 200, 400])
betas = array([5, 10, 15])
temperatures = 1./betas
colors = [plt.cm.jet(i/max(float(len(betas)-1),1)) for i in range(len(betas))]
linestyles = ['-', '--', '-.', ':']
linestyles = ['-', '--', ':']

set_mpl()
fig = plt.figure()
ax_exists = False
ax_nr = 1
for i, a_beta, u in zip(range(len(a_u_beta)), a_u_beta, us):
    if not i in u_inds: continue
    if not ax_exists:
        ax = fig.add_subplot(1,len(u_inds),ax_nr)
        ax.set_ylabel('$A(\omega)$')
        ax_exists = True
    else:
        ax = fig.add_subplot(1,len(u_inds),ax_nr)#,sharey = ax)
        ax.set_yticks(yticks)
        ax.set_yticklabels([])
        ax.set_ylim(yticks.min(), yticks.max())
    ax_nr += 1
    ax.set_xlabel('$\omega$')
    ax.set_title('$U = '+str(u)+'$')
    for a, beta, color, linestyle in zip(a_beta[::-1], betas[::-1], colors, linestyles):
        ax.plot(a[0,:], a[1,:], label = str(1./beta)[:6], color = color, linestyle = linestyle)
    #ax.set_xlim(-.5,2.5)
    ax.set_xlim(-4,4)
    yticks = ax.get_yticks()

plt.legend(loc = 'upper right', title = '$T$', bbox_to_anchor = (.5,1))
plt.tight_layout()
plt.savefig('spectral_function.pdf', dpi=300)
plt.close()
