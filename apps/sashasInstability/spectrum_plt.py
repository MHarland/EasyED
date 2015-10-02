from numpy import load, array, linspace
from matplotlib import pyplot as plt
#from matplotlib_to_latex import set_poster_parameters as set_mpl

a_u_beta = load('spectrum.npy')
mus = [.26,.27,.28]
u_inds = range(3)
betas = array([10, 17, 20])
beta_inds = range(3)

linestyles = ['solid', '--', ':']
temperatures = 1./betas
colors = [plt.cm.jet(i/float(max(len(beta_inds),2)-1)) for i in range(len(beta_inds))]
#set_mpl()
fig = plt.figure()
ax_exists = False
ax_nr = 1
for i, a_beta, mu in zip(range(len(a_u_beta)), a_u_beta, mus):
    if not i in u_inds: continue
    if not ax_exists:
        ax = fig.add_subplot(len(u_inds),1,ax_nr)
        ax_exists = True
    else:
        ax = fig.add_subplot(len(u_inds),1,ax_nr,sharey = ax)
    ax_nr += 1
    ax.set_xlabel('$\omega$')
    ax.set_title('$\mu = '+str(mu)+'$')
    ax.set_ylabel('$\mathrm{Degeneracy}$')
    for a, beta_ind, color, linestyle in zip(a_beta[::-1], beta_inds[::-1], colors, linestyles):
        beta = betas[beta_ind]
        i = 0
        for e, d in zip(*a):
            if i == 0:
                ax.plot([e]*2, [0,d], label = str(1./beta)[:7], color = color, linestyle = linestyle)
            else:
                ax.plot([e]*2, [0,d], color = color, linestyle = linestyle)
            i += 1
    ax.set_xlim(-.02,.05)
    ax.set_ylim(0,5)

plt.legend(loc = 'upper right', title = '$T$')
plt.tight_layout()
plt.savefig('spectrum.pdf', dpi=300)
plt.close()
