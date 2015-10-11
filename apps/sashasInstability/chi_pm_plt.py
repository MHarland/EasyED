from numpy import load, array
from matplotlib import pyplot as plt
#from matplotlib_to_latex import set_poster_parameters as set_mpl

chi_mu_beta = load('chi_pm.npy')
mus = [.26,.27,.28]
inds = [0,1,2]
betas = array([20])
temperatures = 1./betas
colors = [plt.cm.jet(i/max(1,float(len(betas)-1))) for i in range(len(betas))]
linestyles = ['-', '--', ':']
linestyles = ['-']*3
xlims = [0,4]
ylims1 = [0,1.5]
ylims2 = [-.5,.2]
xticklabels = [0,'',1,'',2,'',3,'',4]

#set_mpl()
fig = plt.figure()
ax_exists = False
ax_nr = 1
for i, chi_beta, mu in zip(range(len(chi_mu_beta)), chi_mu_beta, mus):
    if not i in inds: continue
    if not ax_exists:
        ax = fig.add_subplot(1,len(inds),ax_nr)
        ax.set_ylabel('$\chi^{loc}_{+-}(\omega)$')
        ax_exists = True
    else:
        ax = fig.add_subplot(1,len(inds),ax_nr)#, sharey = ax)
        #ax.set_yticks(yticks)
        #ax.set_yticklabels([])
    ax_nr += 1
    ax.set_title('$\mu = '+str(mu)+'$')
    ax.set_xlabel('$\omega$')
    for chi, beta, color, ls in zip(chi_beta[::-1], betas[::-1], colors, linestyles):
        ax.plot(chi[0,:], chi[2,:], label = str(1./beta)[:6], color = color, linestyle = ls)
        #ax.scatter(chi[0, chi[2,:].argmin()], chi[2,:].min(), color = color, marker = '+')
    ax.set_xlim(*xlims)
    #ax.set_ylim(*ylims1)
    yticks = ax.get_yticks()
    #ax.set_xticklabels(xticklabels)
plt.legend(loc = 'upper right', title = '$T$')
plt.tight_layout()
plt.savefig('chi_pm_loc.pdf', dpi=300)
plt.close()

fig = plt.figure()
ax_exists = False
ax_nr = 1
axes = list()
for i, chi_beta, mu in zip(range(len(chi_mu_beta)), chi_mu_beta, mus):
    if not i in inds: continue
    if not ax_exists:
        ax = fig.add_subplot(1,len(inds),ax_nr)
        axes.append(ax)
        ax.set_ylabel('$\chi^{nn}_{+-}(\omega)$')
        ax_exists = True
    else:
        ax = fig.add_subplot(1,len(inds),ax_nr)#,sharey = ax)
        #ax.set_yticks(yticks)
        #ax.set_yticklabels([])
    ax_nr += 1
    ax.set_title('$\mu = '+str(mu)+'$')
    ax.set_xlabel('$\omega$')
    for chi, beta, color, ls in zip(chi_beta[::-1], betas[::-1], colors, linestyles):
        ax.plot(chi[0,:], chi[4,:], label = str(1./beta)[:6], color = color, linestyle = ls)
        #ax.scatter(chi[0, chi[4,:].argmax()], chi[4,:].max(), color = color, marker = '+')
    ax.set_xlim(*xlims)
    #ax.set_ylim(*ylims2)
    yticks = ax.get_yticks()
    #ax.set_xticklabels(xticklabels)
axes[0].legend(loc = 'lower right', title = '$T$')
plt.tight_layout()
plt.savefig('chi_pm_nn.pdf', dpi=300)
plt.close()
