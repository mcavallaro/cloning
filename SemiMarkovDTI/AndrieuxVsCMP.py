import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as spo
from matplotlib import rc
rc('font',**{'family':'serif','serif':['cm']})
rc('text', usetex=True)

def legendre(X,S):
    if len(X) != len(S):
        print 'X and S must have the same len'
        return 0
    lt = []
    J = r_[-10:10:0.01]
    for j in J:
        tmp = max (  [X[i]-S[i]*j for i in range(0,len(S)) ])
        lt.append(tmp)
    return lt

def Easep(s,a,b,g,d):
    return (a+b+g+d-np.sqrt((a+b+g+d)*(a+b+g+d)+4*( (b+g*np.exp(-s))*(a*np.exp(s)+d)-(a+d)*(b+g))))/2


fig = plt.figure()
ax = fig.add_subplot(111)


simul = np.loadtxt('TA.txt')


smin=-1;smax=1.1

fontsize1 = 18
fontsize2 = 14


p21L=0.6; p21R=0.4; p12R=0.5; p12L=0.5

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

k2=1.; k1=1. ; a2=1. ; a1=1.
eNR = [-spo.newton(det, 0.1, args=(s,) ) for s in x]
plt.plot(x,eNR,'kx',label = 'Andrieux-Gaspard solution',markersize=6)#NR k1=k2=a1=a2=1
plt.plot(simul[:,0],-simul[:,1], 'r.', markersize=6,label='cloning')

#k2=1. ; k1=0.01 ; a2=1 ; a1=0.1
#eNR = [-spo.newton(det, 0.1, derdet , args=(s,) ) for s in x]
#plot(x,eNR,'r--', label = 'NR k2=1, k1=0.01, a1=0.1, a2=1')


#k2=1. ; k1=0.01 ; a2=0.1 ; a1=1.
#eNR = [-spo.newton(det, 0.1, derdet , args=(s,) ) for s in x]
#plot(x,eNR,'r.', label = 'NR k2=1, k1=0.01, a1=1, a2=0.1')



x = np.r_[smin:smax:0.0001]



#closed-form solution of Andrieux-Gaspard per k1=k2=0
def e(s):
    return k2*(1-pow(f(s),1/(a2+a1) ))

k2=1.; k1=1. ; a2=1. ; a1=1.
y = [e(s) for s in x]
#plot(x,y,'b-',label='closed form solution k1=k2=1')


#s-Hamiltonial of the corresponding Markov Chain and lowst e-value.

def sH11(s, k1, k2, a1, a2):
    return  - k1 / a1 *(p12L + p12R)

def sH12(s, k1, k2, a1, a2):
    return k2 / a2 *(p21L* np.exp(-s) + p21R)

def sH21(s, k1, k2, a1, a2):
    return k1 / a1 *(p12L * np.exp(s)+ p12R)

def sH22(s, k1, k2, a1, a2):
    return - k2 / a2 *(p21L + p21R)

def em(s, k1, k2, a1, a2):
    matrix = -np.array([[sH11(s, k1, k2, a1, a2) , sH12(s, k1, k2, a1, a2) ],[ sH21(s, k1, k2, a1, a2) , sH22(s, k1, k2, a1, a2) ]])
    return  np.linalg.eig(matrix)[0]


k2=1.; k1=1. ; a2=1. ; a1=1.
Y = [min(em(s, k1, k2, a1, a2)) for s in x]
# plot(x,Y,'r-', label =  'Corresponding Markov Process')

#k2=1. ; k1=0.01 ; a2=1 ; a1=0.1
#Y = [em(s, k1, k2, a1, a2)[0] for s in x]
#plot(x,Y,'g--', label =  'cMP k2=1, k1=0.01, a1=0.1, a2=1')

#k2=1. ; k1=0.01 ; a2=0.1 ; a1=1.
#Y = [min(em(s, k1, k2, a1, a2)) for s in x]
#plot(x,Y,'g-.', label =  'cMP k2=1, k1=0.01, a1=1, a2=0.1')


#effective ASEP
a= k1/a1*p12L
b= k2/a2*p21R
g= k2/a2*p21L
d= k1/a1*p12R
X = np.r_[smin:smax:0.03]
EASEP = [Easep(s, a , b , g, d) for s in X]
#plot(X,EASEP,'o',color='orange', alpha=0.5, label='Effective Exclusion Pro')


print 'The rates of the corresponding Exclusion Process are:'
print 'a=',a,'b=',b,'g=',g,'d=',d

t = ax.text(-0.8, 0.02, r'(b)', ha="center", va="center", size= 20)


plt.axhline(0,ls=':',c='k',linewidth='0.05')
plt.axvline(0,ls=':',c='k',linewidth='0.05')

leg=plt.legend(loc='best', fancybox=True,  frameon=True, prop={'size':fontsize2})
leg.get_frame().set_linewidth(0.0)

#yscale('log')
#ylim([0.000001,0.01])

plt.xlabel(r'$s$')
plt.ylabel(r'$e(s)$')


for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(fontsize2)

plt.xlabel(r'$s$',fontsize=fontsize1)
plt.ylabel(r'$e(s)$',fontsize=fontsize1)


fig.set_size_inches(4*1.2,3*1.2)
plt.tight_layout()

plt.savefig('figure1b.pdf')
