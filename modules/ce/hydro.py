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
mdot=nml2['rhd_quantities']['mdot']
cwd=os.getcwd()
os.chdir(cwd+'/out')
nblk=[]
k=int(sys.argv[1] if len(sys.argv)>=2 else 1)
filenumber=str(k).zfill(5)
filename=filehead+filenumber+'.h5'
print(filename)
t_domain=[]
l_history=[]
while os.path.isfile(filename):
    time,nblocks=read_attr(filename)
    t_domain.append(time)
    x_axis,level_x_axis=block_construct_x_axis(filename,nblocks)
    x_axis_rsun=x_axis/rsun
    dx_x_axis=block_construct_dx(filename,nblocks)
    x_axis_surface=block_construct_xsurface_variable(filename,nblocks,'mesh_x')
    rho=block_construct_cell_variable(filename,nblocks,'rho')
    v=block_construct_cell_variable(filename,nblocks,'vx')
    p=block_construct_cell_variable(filename,nblocks,'pres')
    temp=block_construct_cell_variable(filename,nblocks,'temp')
    nx=np.shape(x_axis)[0]
    fig,axes=plt.subplots(2,figsize=(12,8))
    axes[0].plot(x_axis_rsun,temp,'r-',linewidth=1,label=r'$T_{gas}$')
    axes[0].set_xscale('log')
    axes[0].set_xticks([10,20,100,200,1000,2000,4000])
    axes[0].set_xlim(xmin_rsun,xmax_rsun)
    axes[0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[0].set_xlabel(r'$r/R_{\odot}$')
    axes[0].set_yscale('log')
    axes[0].set_yticks([100,1000,2000,10000,20000])
    axes[0].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[0].set_ylabel('K')
    axes[0].legend(loc=0)


    ax2=axes[1].twinx()
    axes[1].plot(x_axis_rsun,v/1e5,'r-',linewidth=1,label='v')
    axes[1].set_xlim(xmin_rsun,xmax_rsun)
    axes[1].legend(loc=0)
    axes[1].set_xscale('log')
    axes[1].set_xticks([10,20,100,200,1000,2000,4000])
    axes[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[1].set_xlabel(r'$r/R_{\odot}$')
    axes[1].set_ylabel('km/s')
    ax2.plot(x_axis_rsun,rho,'k-',linewidth=1,label='density')
    ax2.legend(loc=5)
    ax2.set_yscale('log')
    ax2.set_ylabel(r'$g\cdot cm^{-3}$')



    plt.tight_layout()
    plt.savefig(cwd+'/explosion'+filenumber+'.png',dpi=100)
    #plt.savefig(cwd+'/explosion'+filenumber+'.pdf',dpi=300)
    plt.close()
    k=k+1
    filenumber=str(k).zfill(5)
    filename=filehead+filenumber+'.h5'
    print(filename)
