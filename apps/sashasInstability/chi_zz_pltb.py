from numpy import load, array, pi
from matplotlib import pyplot as plt

chi_mu_beta = load('chi_zz.npy')
mus = [.26, .27, .28]
betas = array([20])
temperatures = 1./betas
inds = [0,1,2]
colors = [plt.cm.jet(i/2.) for i in range(3)]
linestyles = ['-', '--', ':']
linestyles = ['-']*3
xlims = [0,4]
ylims = [-.6,.6]
xticklabels = [0,'',1,'',2,'',3,'',4]
assert len(betas) == 1, 'Use chi_zz_plt.py for several betas.'

fig = plt.figure()
for j, transition, color  in zip(range(3), ['loc', 'nn', 'nnn'], colors):
    ax_exists = False
    ax_nr = 1
    for i, chi_beta, mu in zip(range(len(chi_mu_beta)), chi_mu_beta, mus):
        if not i in inds: continue
        if not ax_exists:
            ax = fig.add_subplot(1,len(inds),ax_nr)
            ax.set_ylabel('$-4\,\mathrm{Im}\chi_{zz}(\omega)\,/\pi$')
            ax_exists = True
        else:
            ax = fig.add_subplot(1,len(inds),ax_nr)
            ax.set_yticklabels([])
        #ax.set_xticklabels(xticklabels)
        ax.set_ylim(ylims)
        ax_nr += 1
        ax.set_title('$\mu = '+str(mu)+',\, \\beta = '+str(betas[0])+'$')
        ax.set_xlabel('$\omega$')
        for chi in chi_beta[::-1]:
            ax.plot(chi[0,:], -4*chi[1+2*j,:]/pi, label = transition, color = color)
        yticks = ax.get_yticks()
plt.legend(loc = 'upper right')
plt.tight_layout()
plt.savefig('chi_zz_allTrans.pdf', dpi=300)
plt.close()
