import sys
#sys.path.insert(1,'../../scripts')
import numpy as np
import h5py
import f90nml
au=1.496e13
parsec=3.086e+18
G=6.67259e-8
lsun=3.99e33
msun=1.99e33
rsun=6.96e10
rjupiter=7.140e9
mjupiter=1.899e30
yr=3.154e7
amu=1.660539040e-24
arad=7.5646e-15
kb=1.380658e-16
h_planck=6.6260755e-27
c_light=2.9979e10
mearth=5.976e27
day=86400
sigma_sb=5.67051e-5
NA=6.0221367e23
def read_attr(filename):
    f=h5py.File(filename,'r')
    time=f.attrs['time'][0]
    nblocks=f.attrs['nblocks'][0]
    return time,nblocks
def end_point_to_middle_point(x_end,n_end):
    x_middle=np.zeros(n_end-1)
    for i in range(n_end-1):
        x_middle[i]=(x_end[i]+x_end[i+1])/2
    return x_middle
def block_construct_x_axis(filename,nblock):
    f=h5py.File(filename,'r')
    x_axis=np.zeros(0)
    level_x_axis=np.zeros(0)
    for i in range(1,nblock+1):
        groupname='blk'+str(i).zfill(5)
        x_end=f[groupname]['mesh_x']
        n_end=np.shape(x_end)[0]
        block_x=end_point_to_middle_point(x_end,n_end)
        block_level=np.ones(n_end-1)*x_end.attrs['level']
        x_axis=np.concatenate((x_axis,block_x))
        level_x_axis=np.concatenate((level_x_axis,block_level))
    return x_axis,level_x_axis
def block_construct_dx(filename,nblock):
    f=h5py.File(filename,'r')
    dx_x_axis=np.zeros(0)
    for i in range(1,nblock+1):
        groupname='blk'+str(i).zfill(5)
        x_end=f[groupname]['mesh_x']
        n_end=np.shape(x_end)[0]
        block_dx=np.ones(n_end-1)*(x_end[1]-x_end[0])
        dx_x_axis=np.concatenate((dx_x_axis,block_dx))
    return dx_x_axis
def block_construct_cell_variable(filename,nblock,var):
    f=h5py.File(filename,'r')
    v=np.zeros(0)
    for i in range(1,nblock+1):
        groupname='blk'+str(i).zfill(5)
        v_block=f[groupname][var]
        v=np.concatenate((v,v_block[0]))
    return v
def block_construct_cell_variable_multiple(filename,nblock,var,nvar):
    f=h5py.File(filename,'r')
    for i in range(1,nblock+1):
        groupname='blk'+str(i).zfill(5)
        v_block=f[groupname][var][0]
        if (i==1):
            v=v_block
        else:
            v=np.concatenate((v,v_block),axis=0)
    return np.transpose(v)
def block_construct_xsurface_variable(filename,nblock,var):
    f=h5py.File(filename,'r')
    v=np.zeros(0)
    for i in range(1,nblock+1):
        groupname='blk'+str(i).zfill(5)
        v_block=f[groupname][var]
        if (var=='mesh_x'):
            if (i==1):
                v=np.concatenate((v,v_block))
            else:
                nx=np.shape(v_block)[0]
                v=np.concatenate((v,v_block[1:nx]))
        else:
            if (i==1):
                v=np.concatenate((v,v_block[0]))
            else:
                nx=np.shape(v_block)[1]
                v=np.concatenate((v,v_block[0][1:nx]))
    return v
