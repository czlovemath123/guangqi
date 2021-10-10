import sys
sys.path.insert(1,'../../scripts')
from assemble_1d_data import *
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path
plt.rcParams.update({'font.size': 18})
np.set_printoptions(threshold=sys.maxsize)

def read_history(filename):
    nl=len(open(filename).readlines(  ))
    f=open(filename,'r')
    line=f.readline()
    a=[float(s) for s in line.split()]
    da=len(a)
    history=np.zeros((nl,da))
    history[0]=a
    for i in range(nl-1):
        line=f.readline()
        a=[float(s) for s in line.split()]
        history[i+1]=a
    return history,nl



cwd=os.getcwd()
h1,n1=read_history('history1.data')
h1=np.transpose(h1)
h2,n2=read_history('history2.data')
h2=np.transpose(h2)
h3,n3=read_history('history3.data')
h3=np.transpose(h3)
#os.chdir(cwd+'/problem')
nml=f90nml.read('global.data')
os.chdir(cwd)
#lrn_phase=[]
#lrn_temp=[]
#lrn_l=[]
#lrn_r=[]
#f=open('TLR_AT2019zhd_FINAL.dat')
#for i in range(64):
#    line=f.readline()
#    line=[float(s) for s in line.split()]
#    lrn_phase.append(line[0])
#    lrn_temp.append(line[1])
#    lrn_l.append(line[3]*1e37)
#    lrn_r.append(line[5]*214)
#lrn_phase=np.asarray(lrn_phase)
#lrn_temp=np.asarray(lrn_temp)
#lrn_l=np.asarray(lrn_l)
#lrn_r=np.asarray(lrn_r)
#f.close()
h1[10]=h1[10]/day
h2[10]=h2[10]/day
h3[10]=h3[10]/day


fig,ax1=plt.subplots(3,1,figsize=(8,8),sharex=True,squeeze=True)
ax1[0].set_ylabel(r'T [K]')
ax1[0].set_ylim(0,6000)
#ax1[0].plot(h1[10],h1[11],'r-',linewidth=4,label=r'$8.63\times10^{43}\rm{erg}$')
#ax1[0].plot(h2[10],h2[11],'g-',linewidth=2,label=r'$1.71\times10^{44}\rm{erg}$')
#ax1[0].plot(h3[10],h3[11],'b-',linewidth=1,label=r'$4.14\times10^{44}\rm{erg}$')
ax1[0].plot(h1[10],h1[11],'r-',linewidth=4,label=r'$1.72\times10^{44}\rm{erg}$')
ax1[0].plot(h2[10],h2[11],'g-',linewidth=2,label=r'$3.43\times10^{44}\rm{erg}$')
ax1[0].plot(h3[10],h3[11],'b-',linewidth=1,label=r'$8.28\times10^{44}\rm{erg}$')
ax1[0].legend()
#ax1[1].plot(h1[10],h1[12],'r-',linewidth=4,label=r'$2.74\times10^{-3}\ M_{\odot}$')
#ax1[1].plot(h2[10],h2[12],'g-',linewidth=2,label=r'$5.48\times10^{-2}\ M_{\odot}$')
#ax1[1].plot(h3[10],h3[12],'b-',linewidth=1,label=r'$1.37\times10^{-2}\ M_{\odot}$')
ax1[1].plot(h1[10],h1[12],'r-',linewidth=4,label=r'$5.48\times10^{-3}\ M_{\odot}$')
ax1[1].plot(h2[10],h2[12],'g-',linewidth=2,label=r'$1.09\times10^{-2}\ M_{\odot}$')
ax1[1].plot(h3[10],h3[12],'b-',linewidth=1,label=r'$2.74\times10^{-2}\ M_{\odot}$')
ax1[1].legend()
#ax1[1].plot(lrn_phase,lrn_r,'y*',markersize=4)
ax1[1].set_ylabel(r'$R\ [R_{\odot}]$')
ax1[1].set_ylim(5e0,4e3)
ax1[2].plot(h1[10],h1[0],'r-',linewidth=4)
ax1[2].plot(h2[10],h2[0],'g-',linewidth=2)
ax1[2].plot(h3[10],h3[0],'b-',linewidth=1)
#ax1[2].plot(lrn_phase,lrn_l,'y*',markersize=4)
lmax=max(np.max(h1[0]),np.max(h2[0]),np.max(h3[0]))
ax1[2].set_xlim(h1[10][0],h1[10][-1])
ax1[2].set_xlabel('t/day')
ax1[2].set_ylabel(r'$L\ [erg\cdot s^{-1}]$')
ax1[2].set_yscale('log')
ax1[2].set_ylim(2e35,lmax*1.1)
#plt.legend()
#plt.show()
plt.tight_layout()
plt.subplots_adjust(hspace=.0)
plt.savefig(cwd+'/combined.pdf',bbox_inches='tight',dpi=300)
#plt.close()
#plt.show()



