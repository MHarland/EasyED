from numpy import load, array
from matplotlib import pyplot as plt

chi_u_beta = load('chi_pm.npy')
us = [.75,1,3,8]
u_inds = range(3)
betas = array([10, 100, 200, 400])
temperatures = 1./betas
colors = [plt.cm.jet(i/float(len(betas)-1)) for i in range(len(betas))]
linestyles = ['-', '--', '-.', ':']

fig = plt.figure()
ax_exists = False
ax_nr = 1
for i, chi_beta, u in zip(range(len(chi_u_beta)), chi_u_beta, us):
    if not i in u_inds: continue
    if not ax_exists:
        ax = fig.add_subplot(1,len(u_inds),ax_nr)
        ax.set_ylabel('$\mathrm{Im}\,\chi^{loc}_{+-}(\omega)$')
        ax_exists = True
    else:
        ax = fig.add_subplot(1,len(u_inds),ax_nr)#, sharey = ax)
        ax.set_yticks(yticks)
        ax.set_ylim(yticks.min(), yticks.max())
        ax.set_yticklabels([])
    ax_nr += 1
    ax.set_title('$U = '+str(u)+'$')
    ax.set_xlabel('$\omega$')
    for chi, beta, color, ls in zip(chi_beta[::-1], betas[::-1], colors, linestyles):
        ax.plot(chi[0,:], chi[1,:], label = str(1./beta)[:6], color = color, linestyle = ls)
        ax.scatter(chi[0, chi[1,:].argmin()], chi[1,:].min(), color = color, marker = '+')
    #ax.set_xlim(-.5,2.5)
    yticks = ax.get_yticks()
    ax.set_ylim(-3,0)
plt.legend(loc = 'upper right', title = '$T$')
plt.tight_layout()
plt.savefig('chi_pm_loc.pdf', dpi=300)
plt.close()

fig = plt.figure()
ax_exists = False
ax_nr = 1
for i, chi_beta, u in zip(range(len(chi_u_beta)), chi_u_beta, us):
    if not i in u_inds: continue
    if not ax_exists:
        ax = fig.add_subplot(1,len(u_inds),ax_nr)
        ax.set_ylabel('$\mathrm{Im}\,\chi^{nn}_{+-}(\omega)$')
        ax_exists = True
    else:
        ax = fig.add_subplot(1,len(u_inds),ax_nr)#,sharey = ax)
        ax.set_yticks(yticks)
        ax.set_yticklabels([])
        ax.set_ylim(yticks.min(), yticks.max())
    ax_nr += 1
    ax.set_title('$U = '+str(u)+'$')
    ax.set_xlabel('$\omega$')
    for chi, beta, color, ls in zip(chi_beta[::-1], betas[::-1], colors, linestyles):
        ax.plot(chi[0,:], chi[3,:], label = str(1./beta)[:6], color = color, linestyle = ls)
        ax.scatter(chi[0, chi[3,:].argmax()], chi[3,:].max(), color = color, marker = '+')
    #ax.set_xlim(-.5,2.5)
    yticks = ax.get_yticks()
plt.legend(loc = 'lower right', title = '$T$')
plt.tight_layout()
plt.savefig('chi_pm_nn.pdf', dpi=300)
plt.close()
