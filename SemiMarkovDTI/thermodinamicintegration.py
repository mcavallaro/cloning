from numpy import *
import  math, matplotlib.pyplot as plt 
from matplotlib import rc
import sys
rc('font',**{'family':'serif','serif':['cm']})
rc('text', usetex=True)


if len(sys.argv) != 3:
    print 'Usage: python Thermodynamicintegration.py [file_input] [file_output.pdf] ')
    sys.exit(1)


def SCGF(s,a,b,g,d):
    SQRT = sqrt((a+b+g+d)*(a+b+g+d)+4*((b+g*exp(s))*(a*exp(-s)+d)-(a+d)*(b+g)))
    return (a+b+g+d -SQRT)/2


data = loadtxt(sys.argv[0])
output_file=sys.argv[1]

al=0.1
be=0.2
ga=0
de=0
deltas = 0.1


sv = data[:,0]
eigen = data[:,1]
meanobsv = data[:,3]

ll = len(sv)

fig = plt.figure()
ax = fig.add_subplot(111)

length = len(meanobsv)
observable = [0 for i in range(length)]
theoretical = [  SCGF(s,al,be,ga,de)   for s in sv ]

zeropoint=int(length/2)

#integrate
for k in range(1,zeropoint+1):
	observable[zeropoint+k] = observable[zeropoint+k-1] + meanobsv[zeropoint+k]*deltas
	observable[zeropoint-k] = observable[zeropoint-k+1] - meanobsv[zeropoint-k]*deltas

plt.axhline(0,ls=':',c='k',linewidth=0.05)
plt.axvline(0,ls=':',c='k',linewidth=0.05)

plt.plot(sv, -eigen,'y.',label="Direct evaluation")
plt.plot(sv, observable,'g.', label="Thermodynamic integration")
plt.plot(sv, theoretical, 'k+',label='theoretical')

plt.legend(loc='best')

plt.xlabel('$s$',fontsize=20)
plt.ylabel('$e(s)$',fontsize=20)

for label in ax.get_xticklabels() + ax.get_yticklabels():
	label.set_fontsize(15)

fig = plt.gcf()
fig.set_size_inches(4*1.5,3*1.5)
plt.tight_layout()

plt.savefig(output_file)
plt.close()

