from numpy import load, array, pi
from matplotlib import pyplot as plt
#from matplotlib_to_latex import set_poster_parameters as set_mpl

chi_u_beta = load('chi_c.npy')
betas = array([10,15,20,25,30])
us = [2.75,3,3.25]
temperatures = 1./betas
inds = [0,1,2]
colors = [plt.cm.jet(i/float(max(len(betas), 2)-1)) for i in range(len(betas))]
xlims = [0,8]
ylims = [-.1,.4]
xticklabels = [0,'',1,'',2,'',3,'',4]

#set_mpl()
fig = plt.figure()
for j, transition, ls in zip(range(2), ['loc', 'nn'], ['solid', 'dashed']):
    ax_exists = False
    ax_nr = 1
    for i, chi_beta, u in zip(range(len(chi_u_beta)), chi_u_beta, us):
        if not i in inds: continue
        if not ax_exists:
            ax = fig.add_subplot(1,len(inds),ax_nr)
            ax.set_ylabel('$-\mathrm{Im}\chi_{c}(\omega)\,/\pi$')
            ax_exists = True
        else:
            ax = fig.add_subplot(1,len(inds),ax_nr)#, sharey = ax)
            ax.set_yticklabels([])
        ax.set_ylim(ylims)
        ax_nr += 1
        ax.set_title('$U = '+str(u)+'$')
        ax.set_xlabel('$\omega$')
        for chi, color in zip(chi_beta[::-1], colors):
            ax.plot(chi[0,:], -chi[1+2*j,:]/pi, label = transition, color = color, ls = ls)
        yticks = ax.get_yticks()
#plt.legend(loc = 'upper right')
plt.tight_layout()
plt.savefig('chi_c_allTransf.pdf', dpi=300)
plt.close()
