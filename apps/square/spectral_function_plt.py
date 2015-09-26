from numpy import load, array
from matplotlib import pyplot as plt

a_u_beta = load('spectral_function.npy')
us = [1,8,12]
betas = array([1,5,10])
u_inds = range(3) # indices of us that will be plotted
temperatures = 1./betas
colors = [plt.cm.jet(i/float(len(betas)-1)) for i in range(len(betas))]
linestyles = ['-', '--', ':']

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
    #ax.set_xlim(-.2,1.2)
    #ax.set_ylim(0,3)
    yticks = ax.get_yticks()

plt.legend(loc = 'upper right', title = '$T$')
plt.tight_layout()
plt.savefig('spectral_function.pdf', dpi=300)
plt.close()
