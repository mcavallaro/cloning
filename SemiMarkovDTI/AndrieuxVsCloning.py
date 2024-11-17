import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as spo
from matplotlib import rc
rc('font',**{'family':'serif','serif':['cm']})
rc('text', usetex=True)


fig = plt.figure()
ax = fig.add_subplot(111)

# cloning1 = np.loadtxt('SCGF_1.txt')
cloning2 = np.loadtxt('SCGF_2.txt')
# cloning3 = np.loadtxt('SCGF_3.txt')


smin=-1;smax=1.1

fontsize1 = 18
fontsize2 = 14

p21L=0.6; p21R=0.4; p12R=0.5; p12L=0.5 #(alpha, beta, gamma, delta)

x = np.r_[smin:smax:0.1]

def f(s):
    return (p21L * np.exp(-s) + p21R)*(p12L *np.exp(s) + p12R)

A  = 0.405

#Newton-Raphson solution of Andrieux-Gaspard per k1=k2=0
# the SCGF is the value of c such that det(c)==0
# check if C is conjugated to the time
def det(c,s):
    return pow((1+c/k1),a1)*pow((1+c/k2),a2)  - f(s)

def derdet(c,s):
    return ((a1 * pow((1+c/k1),(-1+a1)) * pow((1+c/k2),a2)) )/k1 + (a2*pow((1+c/k1),a1)*pow( (1+c/k2),(-1+a2)) ) /k2

# k2=1.; k1=1. ; a2=1. ; a1=1.
# eNR1 = [-spo.newton(det, 0.1, args=(s,) ) for s in x]
# plt.plot(x, eNR1, 'kx',label = 'Andrieux-Gaspard solution',markersize=6) #NR k1=k2=a1=a2=1
# plt.plot(cloning1[:,0],-cloning1[:,1], 'r.', markersize=6,label='cloning')

k2=1. ; k1=0.01 ; a2=1. ; a1=0.1
eNR2 = [-spo.newton(det, 0 , derdet, args=(s,),maxiter=50 ) for s in x]
#plt.plot(x,eNR2,'kx', label = 'Solution of Ref.~[38]')
plt.plot(x,eNR2,'kx', label = 'Sol. of Andrieux and Gaspard (2008)')
plt.plot(cloning2[:,0],-cloning2[:,1], 'r.', markersize=6, label='Cloning')


# k2=1. ; k1=0.01 ; a2=0.1 ; a1=1.
# eNR3 = [-spo.newton(det, 0.1, derdet , args=(s,) ) for s in x]
# plt.plot(x,eNR3,'kx', label = 'NR k2=1, k1=0.01, a1=1, a2=0.1')
# plt.plot(cloning3[:,0],-cloning3[:,1], 'r.', markersize=6,label='cloning k2=1, k1=0.01, a1=0.1, a2=1')



t = ax.text(0.05, 0.85,  '(a)', size= 20, transform=ax.transAxes, backgroundcolor='white')


plt.axhline(0,ls=':',c='k',linewidth=0.05)
plt.axvline(0,ls=':',c='k',linewidth=0.05)

leg=plt.legend(loc='best', fancybox=True,  frameon=True, prop={'size':fontsize2-0.5})
leg.get_frame().set_linewidth(0.0)

#yscale('log')
#ylim([0.000001,0.01])

plt.xlabel(r'$s$')
plt.ylabel(r'$e(s)$')


for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(fontsize2)

plt.xlabel(r'$s$',fontsize=fontsize1)
plt.ylabel(r'$e(s)$',fontsize=fontsize1)

X=1.2
#X=5
fig.set_size_inches(4*X,3*X)
plt.tight_layout()

plt.savefig('figure.pdf')
