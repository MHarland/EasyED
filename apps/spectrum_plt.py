from numpy import load, array, linspace
from matplotlib import pyplot as plt
#from matplotlib_to_latex import set_poster_parameters as set_mpl

a_u_beta = load('spectrum.npy')
us = [2.75,3,3.25]
u_inds = [0,1,2]
betas = array([10,15,20,25,30])
#betas = array([10,69,70,86,87])
beta_inds = range(len(betas))
linestyles = ['-','-', '--', '-.', ':']

temperatures = 1./betas
colors = [plt.cm.jet(i/float(max(len(beta_inds),2)-1)) for i in range(len(beta_inds))]
#set_mpl()
fig = plt.figure()
axes = list()
ax_exists = False
ax_nr = 1
for i, a_beta, u in zip(range(len(a_u_beta)), a_u_beta, us):
    if not i in u_inds: continue
    if not ax_exists:
        ax = fig.add_subplot(len(u_inds),1,ax_nr)
        ax_exists = True
    else:
        ax = fig.add_subplot(len(u_inds),1,ax_nr,sharey = ax)
    axes.append(ax)
    if not ax_nr == len(u_inds):
        ax.set_xticklabels([])
    else:
        ax.set_xlabel('$\omega$')
    ax_nr += 1
    ax.text(.05,.95,'$U = '+str(u)+'$', transform = ax.transAxes, ha = 'left', va = 'top')
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
    #ax.set_xlim(-.02,1.5)
    #ax.set_ylim(0,10)

axes[0].legend(loc = 'upper right', title = '$T$')
plt.tight_layout()
plt.savefig('spectrum.pdf', dpi=300)
plt.close()
