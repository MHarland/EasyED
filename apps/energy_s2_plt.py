from numpy import load, array, linspace
from matplotlib import pyplot as plt
#from matplotlib_to_latex import set_poster_parameters as set_mpl

a_u_beta = load('energy_s2.npy')
us = [3]
u_inds = [0]
betas = array([10])
beta_inds = [0]

temperatures = 1./betas
colors = [plt.cm.jet(i/float(max(len(beta_inds),2)-1)) for i in range(len(beta_inds))]
#set_mpl()
fig = plt.figure()
ax_exists = False
ax_nr = 1
for i, a_beta, u in zip(range(len(a_u_beta)), a_u_beta, us):
    if not i in u_inds: continue
    if not ax_exists:
        ax = fig.add_subplot(len(u_inds),1,ax_nr)
        ax_exists = True
    else:
        ax = fig.add_subplot(len(u_inds),1,ax_nr,sharey = ax)
    ax_nr += 1
    ax.set_xlabel('$<S^2_{tot}>$')
    ax.set_title('$U = '+str(u)+'$')
    ax.set_ylabel('$E$')
    for a, beta_ind, color in zip(a_beta[::-1], beta_inds[::-1], colors):
        beta = betas[beta_ind]
        i = 0
        for s2, e, deg in zip(*a):
            if i == 0:
                ax.scatter(s2, e, label = str(1./beta)[:7], color = color, marker = '_', zorder = 2)
            else:
                ax.scatter(s2, e, color = color, marker = '_', zorder = 2)
            ax.text(s2, e, str(int(deg)), color = 'gray')
            i += 1
    s2_heisenberg = [(1./u)*s*(s+1) for s in range(100)]
    for s2 in s2_heisenberg:
        ax.plot([s2]*2,[-.05*a_u_beta[:,:,1].max(),1.05*a_u_beta[:,:,1].max()], color = 'gray', linestyle = 'dashed', zorder = 1)
    ax.set_xlim(-.05*a_u_beta[:,:,0].max(),1.05*a_u_beta[:,:,0].max())
    ax.set_ylim(-.05*a_u_beta[:,:,1].max(),1.05*a_u_beta[:,:,1].max())
    #ax.set_ylim(-.01,1)
    xticks = list()
    for s2 in a[0,:]:
        is_element = False
        for tick in xticks:
            if abs(s2-tick) < 10**(-10):
                is_element = True
        if not is_element:
            xticks.append(s2)
    ax.set_xticks(xticks)

plt.legend(loc = 'upper right', title = '$T$')
plt.tight_layout()
plt.savefig('energy_s2.pdf', dpi=300)
plt.close()
