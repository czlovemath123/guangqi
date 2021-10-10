import sys
sys.path.insert(1,'../../scripts')
from assemble_1d_data import *
from math import *
import matplotlib
import matplotlib.pyplot as plt
import os.path
plt.rcParams.update({'font.size': 20})
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

def blackbody(teff,wavelength):
    e=2*h_planck*c_light*c_light/(pow(wavelength,5)*(exp(h_planck*c_light/wavelength/kb/teff)-1))
    return e

def iband(teff):
    if (teff<1000):
        e=1e-5
    else:
        nps=100
        e=0
        wl0=(806-149)*1e-7
        wl1=(806+149)*1e-7
        dl=(wl1-wl0)/nps
        for i in range(nps):
            wl=wl0+dl*i
            p=blackbody(teff,wl)
            e=e+p*dl
        e=e/sigma_sb/pow(teff,4)
    return e


cwd=os.getcwd()
os.chdir(cwd+'/out')
h,n=read_history('history.data')
h=np.transpose(h)
os.chdir(cwd)
nml=f90nml.read('global.data')
tend=nml['meshinfo']['tfinal']
tend=4e6
t0=np.linspace(0,tend,num=n+1)
t=np.zeros(n)
t=t0[1:n+1]/day

fig,ax1=plt.subplots(figsize=(16,8))
ln1=ax1.plot(t,h[1],'k-',linewidth=1,label=r'$e_{g}$')
ln2=ax1.plot(t,h[2],'r-',linewidth=1,label=r'$e_{k}$')
ln3=ax1.plot(t,h[3],'b-',linewidth=1,label=r'$e_{rad}$')
ax1.set_xlim(t[0],t[n-1])
ax1.set_xlabel('day')
ax1.set_yscale('log')
ax1.set_ylabel(r'$L_{\odot}\times day$')
lns=ln1+ln2+ln3
labs=[l.get_label() for l in lns]
ax1.legend(lns,labs,loc=0)
#plt.show()
plt.savefig(cwd+'/energy_budget.png',dpi=300)
#plt.savefig(cwd+'/evolution.pdf',dpi=300)




fig,ax1=plt.subplots(figsize=(16,8))
ln1=ax1.plot(t,h[0],'k-',linewidth=2,label='L')
ln2=ax1.plot(t,-h[9],'r.',markersize=1,label=r'$-\dot{e}_{total}$')
ln3=ax1.plot(t,-h[10],'b.',markersize=1,label=r'$-\dot{e}_{total}+\dot{e}_{rad}$')
ax1.set_xlim(t[0],t[n-1])
ax1.set_xlabel('day')
ax1.set_ylabel(r'$L/L_{\odot}$')
ax1.set_yscale('log')
ax1.set_ylim(1e0,3e5)
lns=ln1+ln2+ln3
labs=[l.get_label() for l in lns]
ax1.legend(lns,labs,loc=0)
#plt.show()
plt.savefig(cwd+'/dedt.png',dpi=300)



fig,ax1=plt.subplots(figsize=(16,8))
ax1.plot(t,h[11],'k-',linewidth=1,label='Teff')
ax1.set_xlim(t[0],t[n-1])
ax1.set_xlabel('day')
ax1.set_ylabel(r'T [K]')
ax1.set_ylim(300,7000)
plt.legend()
#plt.show()
plt.savefig(cwd+'/teff.png',dpi=300)

#nps=500
#wl=np.zeros(nps)
#specturm=np.zeros(nps)
#wl0=100*1e-7
#dlog=3/nps
#teff=4000
#for i in range(nps):
#    wl[i]=wl0*pow(10,dlog*i)
#    specturm[i]=blackbody(teff,wl[i])
#plt.plot(wl*1e7,specturm)
#plt.yscale('log')
#plt.xscale('log')
#plt.show()



#bi=np.zeros(n)
#for i in range(n):
#    teff=h[11][i]
#    bi[i]=iband(teff)*h[0][i]
#plt.plot(t,bi)
#plt.yscale('log')
#plt.ylim(1e-1,1e5)
##plt.xlim(0,40)
#plt.show()
