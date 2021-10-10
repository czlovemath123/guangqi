import sys
sys.path.insert(1,'../../scripts')
from assemble_1d_data import *
from math import *
import matplotlib
import matplotlib.pyplot as plt
import os.path
plt.rcParams.update({'font.size': 16})
np.set_printoptions(threshold=sys.maxsize)
def total_mass(x,rho,dx,nx):
    m=0
    for i in range(nx):
        m=m+4*pi*rho[i]*x[i]*x[i]*dx
    return m
def find_teff(tau,temp,nx):
    for i in range(nx-1):
        if (tau[i]>1.0/3.0 and tau[i+1]<1.0/3.0):
            break
    teff=temp[i]
    return teff
def find_teff_xindex(tau,nx):
    for i in range(nx-1):
        if (tau[i]>2.0/3.0 and tau[i+1]<2.0/3.0):
            break
    return i
f_var=open('output_var_info.dat','r')
n_var=int(f_var.readline())
print(n_var)
var_list=[]
for i in range(n_var):
    line=f_var.readline()
    var_list.append(line.rstrip())
print(var_list)
line=f_var.readline()
f_var.close()
filehead=line.rstrip()
nml=f90nml.read('global.data')
xmin=nml['meshinfo']['n_domain'][0]
xmax=nml['meshinfo']['n_domain'][1]
xmin_rsun=xmin/rsun
xmax_rsun=xmax/rsun
nframe=nml['meshinfo']['nframe']
tfinal=nml['meshinfo']['tfinal']
dt=tfinal/nframe
nml2=f90nml.read('problem.data')
#mdot=nml2['rhd_quantities']['mdot']
cwd=os.getcwd()
data_dir='/'+sys.argv[1] if len(sys.argv)>=2 else '/out'
fig_dir='/'+sys.argv[2] if len(sys.argv)>=3 else ''
os.chdir(cwd+data_dir)
print(data_dir,fig_dir)
nblk=[]
k=1
filenumber=str(k).zfill(5)
filename=filehead+filenumber+'.h5'
print(filename)
t_domain=[]
while os.path.isfile(filename):
    time,nblocks=read_attr(filename)
    t_domain.append(time)
    x_axis,level_x_axis=block_construct_x_axis(filename,nblocks)
    x_axis_rsun=x_axis/rsun
    dx_x_axis=block_construct_dx(filename,nblocks)
    x_axis_surface=block_construct_xsurface_variable(filename,nblocks,'mesh_x')
    x_axis_surface_rsun=x_axis_surface/rsun
    rho=block_construct_cell_variable(filename,nblocks,'rho')
    v=block_construct_cell_variable(filename,nblocks,'vx')
    p=block_construct_cell_variable(filename,nblocks,'pres')
    Erad=block_construct_cell_variable(filename,nblocks,'Erad')
    Frad=block_construct_xsurface_variable(filename,nblocks,'Fradx')
    kp=block_construct_cell_variable(filename,nblocks,'Planck')
    kr=block_construct_cell_variable(filename,nblocks,'Rosseland')
    temp=block_construct_cell_variable(filename,nblocks,'temp')
    entropy=block_construct_cell_variable(filename,nblocks,'entropy')
    h2=block_construct_cell_variable(filename,nblocks,'H2')
    hi=block_construct_cell_variable(filename,nblocks,'HI')
    hii=block_construct_cell_variable(filename,nblocks,'HII')
    erad_temp=np.power(Erad/arad,0.25)
    nx=np.shape(x_axis)[0]
    tau_planck=1e-5*np.ones(nx)
    tau_rosseland=1e-5*np.ones(nx)
    tau_planck=1e-5*np.ones(nx)
    mfp=np.zeros(nx)
    print(Frad[nx-2]/c_light/Erad[nx-2])
    for i in range(nx-1):
        tau_planck[nx-i-2]=tau_planck[nx-i-1]+(dx_x_axis[nx-i-2]*kp[nx-i-2]*rho[nx-i-2]+dx_x_axis[nx-i-1]*kp[nx-i-1]*rho[nx-i-1])/2
    for i in range(nx-1):
        tau_rosseland[nx-i-2]=tau_rosseland[nx-i-1]+(dx_x_axis[nx-i-2]*kr[nx-i-2]*rho[nx-i-2]+dx_x_axis[nx-i-1]*kr[nx-i-1]*rho[nx-i-1])/2
        tau_planck[nx-i-2]=tau_planck[nx-i-1]+(dx_x_axis[nx-i-2]*kp[nx-i-2]*rho[nx-i-2]+dx_x_axis[nx-i-1]*kp[nx-i-1]*rho[nx-i-1])/2
    for i in range(nx+1):
        Frad[i]=Frad[i]*4*pi*x_axis_surface[i]*x_axis_surface[i]/lsun
    for i in range(nx):
        mfp[i]=1.0/(rho[i]*kr[i])
    mfp=mfp/rsun



    fig,axes=plt.subplots(4,1,figsize=(8,12),sharex=True,squeeze=True)
    ln1=axes[0].plot(x_axis_surface_rsun,Frad,'k-',linewidth=1,label=r'$L$')
    axes[0].set_xlim(xmin_rsun,xmax_rsun)
    axes[0].set_xscale('log')
    axes[0].set_xticks([10,20,100,200,1000,2000])
    axes[0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[0].set_xlabel(r'$r\ [R_{\odot}]$')
    axes[0].set_ylabel(r'$L\ [L_{\odot}]$')
    axes[0].legend(loc=0)


    ax2=axes[1].twinx()
    ln1=axes[1].plot(x_axis_rsun,rho,'k-',linewidth=1,label=r'$\rho$')
    axes[1].set_xlim(xmin_rsun,xmax_rsun)
    axes[1].set_xscale('log')
    axes[1].set_xticks([10,20,100,200,1000,2000])
    axes[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[1].set_xlabel(r'$r\ [R_{\odot}]$')
    axes[1].set_ylabel(r'$\rho\ [g\cdot cm^{-3}]$')
    axes[1].set_yscale('log')
    #ax2.plot(x_axis_rsun,kr,'r-.',linewidth=1,label=r'$\kappa_{Rosseland}$')
    #ln2=ax2.plot(x_axis_rsun,tau_rosseland,'r-',linewidth=1,label=r'$\tau_{Rosseland}$')
    #ln2=ax2.plot(x_axis_rsun,tau_planck,'r-.',linewidth=1,label=r'$\tau_{Planck}$')
    #ln2=ax2.plot(x_axis_rsun,kr,'r-',linewidth=1,label=r'$\kappa_{Planck}$')
    ln2=ax2.plot(x_axis_rsun,mfp/x_axis_rsun,'r-',linewidth=1,label=r'$\lambda_{mfp}/r$')
    ax2.axhspan(1,10,alpha=0.5,color='grey')
    lns=ln1+ln2
    labs=[l.get_label() for l in lns]
    axes[1].legend(lns,labs,loc=0)
    ax2.set_ylabel(r'$\lambda_{mfp}/r$')
    ax2.set_yscale('log')

    ln1=axes[2].plot(x_axis_rsun,erad_temp,'b-',linewidth=1,label=r'$T_{rad}$')
    ln2=axes[2].plot(x_axis_rsun,temp,'r-',linewidth=1,label=r'$T_{gas}$')
    axes[2].set_xscale('log')
    axes[2].set_xticks([10,20,100,200,1000,2000])
    axes[2].set_xlim(xmin_rsun,xmax_rsun)
    axes[2].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[2].set_xlabel(r'$r/R_{\odot}$')
    axes[2].set_yscale('log')
    axes[2].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[2].set_ylabel('T [K]')
    axes[2].set_yscale('log')
    lns=ln1+ln2
    labs=[l.get_label() for l in lns]
    axes[2].legend(lns,labs,loc=0)



    ln1=axes[3].plot(x_axis_rsun,hi,'r-',linewidth=1,label=r'$\chi_{\rm{H}}$')
    ln2=axes[3].plot(x_axis_rsun,hii,'k-',linewidth=1,label=r'$\chi_{\rm{H}^{+}}$')
    ln3=axes[3].plot(x_axis_rsun,h2,'b-',linewidth=1,label=r'$\chi_{\rm{H_{2}}}$')
    axes[3].legend(lns,labs,loc=0)
    axes[3].set_xlim(xmin_rsun,xmax_rsun)
    axes[3].set_xscale('log')
    axes[3].set_xticks([10,20,100,200,1000,2000])
    axes[3].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[3].set_xlabel(r'$r\ [R_{\odot}]$')
    axes[3].set_ylabel(r'$\chi$')
    lns=ln1+ln2+ln3
    labs=[l.get_label() for l in lns]
    axes[3].legend(lns,labs,loc=0)
    plt.tight_layout()
    plt.subplots_adjust(hspace=.0)
    plt.savefig(cwd+fig_dir+'/explosion'+filenumber+'.png',bbox_inches='tight',dpi=100)
    #plt.savefig(cwd+'/explosion'+filenumber+'.pdf',dpi=300)
    plt.close()
    k=k+1
    filenumber=str(k).zfill(5)
    filename=filehead+filenumber+'.h5'
    print(filename)
