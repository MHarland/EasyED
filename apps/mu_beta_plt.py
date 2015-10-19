from numpy import load, array
from matplotlib import pyplot as plt

us = [1,1.5,2,2.75,3,3.25,4,8]
fnames = ['beta_mu_u'+str(u)+'.npy' for u in us]
colors = [plt.cm.jet(float(i)/max(float(len(us)-1),1)) for i in range(len(us))]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for u, fname, color in zip(us, fnames, colors):
    beta_mu = load(fname)
    ax.plot(beta_mu[0,:], beta_mu[1,:], label = str(u), color = color)
ax.set_xlabel('$\\beta$')
ax.set_ylabel('$\\mu$')
#ax.set_xlim([0,50])

plt.legend(title = '$U$')
plt.tight_layout()
plt.savefig('mu_beta.pdf', dpi = 300)
plt.close()
