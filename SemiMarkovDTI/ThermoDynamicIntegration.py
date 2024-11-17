from numpy import *
import  math, matplotlib.pyplot as plt 

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['cm']})
rc('text', usetex=True)

name='N210000_5000_0.1_0.2_0_0_c=-2_2_i=1.txt'
matrix = loadtxt(name)

#name_left='N210000_5000_0.1_0.2_0_0_c=-2_2_i=1.txt'
nameleft='lefttxt'
matrix_left = loadtxt(nameleft)


#name='N210000_5000_0.1_0.2_0_0_c=0.05_-2_2_i=1.txt';
#matrix2 = loadtxt(name)

name='ZRP_right_left'
namefile=name+'.pdf'

al=0.1
be=0.2
c=100000

#la quarta colona e' meanobs
eigen = matrix[:,1]
meanobsv = matrix[:,2]

eigen_left = matrix_left[:,1]
meanobsv_left = matrix_left[:,2]

#eigen2 = matrix2[:,1]
#meanobsv2 = matrix2[:,2]

sv = matrix[:,0]
ll = len(sv)

fig = plt.figure()
ax = fig.add_subplot(111)

deltas = 0.01
#adesso lavoriamo un po con meanobsv
length = len(meanobsv)
observable = [0 for i in range(length)]
observable_left = [0 for i in range(length)]
#observable2 = [0 for i in range(length)]

theoretical = [  -(exp(-s)-1)*al   for s in sv ]
#theoretical_2 = [  -(exp(-s)-1)*(al*be/(be+ga))+(1-exp(s))*(ga*de/(be+ga)) for s in sv ]
#theoretical_3 = [ al + be +ga + de- 2 * sqrt((al* exp(-s) +de)*(be+ga*exp(s))) for s in sv ]

zeropoint=int(length/2)

#the mean value must be divided by the cloningfactor
#for k in range(len(meanobsv)):
#    meanobsv[k]=meanobsv[k]/Clov[k]

#integrate
for k in range(1,zeropoint):
	observable[zeropoint+k] = observable[zeropoint+k-1] + meanobsv[zeropoint+k]*(deltas)
	observable[zeropoint-k] = observable[zeropoint-k+1] - meanobsv[zeropoint-k]*(deltas)
	
for k in range(1,zeropoint):
	observable_left[zeropoint+k] = observable_left[zeropoint+k-1] + meanobsv_left[zeropoint+k]*(deltas)
	observable_left[zeropoint-k] = observable_left[zeropoint-k+1] - meanobsv_left[zeropoint-k]*(deltas)

#for k in range(1,zeropoint):
#    observable2[zeropoint+k] = observable2[zeropoint+k-1] + meanobsv2[zeropoint+k]*(deltas)
#    observable2[zeropoint-k] = observable2[zeropoint-k+1] - meanobsv2[zeropoint-k]*(deltas)
#

plt.axhline(0,ls=':',c='k',linewidth=0.05)
plt.axvline(0,ls=':',c='k',linewidth=0.05)

#eta = sqrt( ((be+ga)*(be+ga) -be*de -al*ga)*((be+ga)*(be+ga) -be*de -al*ga) - 4*al*be*ga*de)
#plt.axvline(log(al/(be+ga-de)),ls=':',c='b',linewidth='0.05')
#plt.axvline(log( ((be+ga)*(be+ga) - al*ga - be*de+eta)/(2*ga*de) ),ls=':',c='b',linewidth='0.05')

ccc = [c for i in sv]



#plt.plot(sv, observable,'r-', label="Thermodynamic integration")
plt.plot(sv, -eigen,'y-',label="Direct simulation")

#plt.plot(sv, observable_left,'g-', label="Thermodynamic integration input curr")
plt.plot(sv, -eigen_left,'b-', label="Direct simulation input curr")

#
#plt.plot(sv, observable2,'-',color="blue", label="Thermodynamic integration ens=5000")
#plt.plot(sv, -eigen2,'g-', label="Direct simulation ens=5000")
#

plt.plot(sv, theoretical, 'k-',label=r'$\lambda(1-e^{-s})$') #linewidth=0.05,

plt.plot(sv, ccc, ls=':', c='k',linewidth=0.05)

plt.legend(loc='best', fancybox=True, prop={'size':5})


textstr1 = r'$\lambda=0.1$'
textstr2 = r'$\mu(n) = 0.2 \cdot n$'
textstr3 = r'$c = 0.05$'
textstr = textstr1 + '\n' +textstr2 + '\n' + textstr3

#plt.axvline(s_crit(al,c),ls=':',c='r', label=r'$S_{crit,c}=log \frac{\alpha }{ (\alpha-c)}$')


plt.xlabel('s',fontsize=20)
plt.ylabel('e(s)',fontsize=20)
plt.title('SCGF, right current, ordinary ZRP, ind. particles')


for label in ax.get_xticklabels() + ax.get_yticklabels():
	label.set_fontsize(15)

# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=1.0)
# place a text box in upper left in axes coords
plt.text(0.25, -0.13, textstr, fontsize=16, verticalalignment='top', bbox=props)

#xtickss=list([-1,-0.5,0,1,1.5]) + list([s_crit(al,c)])
#ytickss=list([-0.2,-0.1,0,0.1,0.2]) + list([c])

str_s1=r'$s_{1}$' #+str(s_2(al,be,c)

#plt.xticks(xtickss, ("-1","-0.5","0","1","1.5",str_s1) )
#plt.yticks(ytickss, ("-0.2","-0.1","0", "0.1", "0.2", r'c') )

plt.ylim(-0.5,0.2)
#plt.ylim(-1,0.3)
#plt.ylim(0,c)

fig = plt.gcf()
fig.set_size_inches(4*1.5,3*1.5)
#plt.tight_layout()

#plt.show()

plt.savefig(namefile)
plt.close()







