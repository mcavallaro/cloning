from pylab import *
from matplotlib import rc
rc('font',**{'family':'serif','serif':['cm']})
rc('text', usetex=True)

fontsize1 = 18
fontsize2 = 14

nome = 'RingDiscr_1000_1000_5_aggr'

simul = loadtxt('directory/'+nome+'.txt')


fig = plt.figure()


ax = fig.add_subplot(111)

X = np.r_[-1:1:0.01]
Y = [ -log(0.6 *exp(-s)+0.4*exp(s))   for s in X]


ax.plot(X,Y, 'r-', label =  r'$-\ln(0.6 e^{-s} + 0.4 e^{s})$')
ax.plot(simul[:,0],-simul[:,1], 'b.', label=r'Discrete-time cloning')

ax.axvline(0,color = 'k',ls=':', lw = 0.2)
ax.axhline(0,color = 'k',ls=':', lw = 0.2)
leg =legend(loc = 'lower center',frameon=True, prop={'size':fontsize2})
leg.get_frame().set_linewidth(0.0)

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(fontsize2)


#ylim([0.000001,0.01])

xlabel(r'$s$',fontsize=fontsize1)
ylabel(r'$e(s)$',fontsize=fontsize1)

fig.set_size_inches(4*1.5,3*1.5)
plt.tight_layout()

savefig(nome+'.pdf')
