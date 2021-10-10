module datastructure
#include <petsc/finclude/petscksp.h>
use mpi
use petscksp
implicit none

real(8), parameter, public ::    &
pi=3.141592653589793238462643383279502884197169399375105d0,  &
e_exp=2.718281828459d0

real(8), parameter, public ::    &
gr=6.67259d-8,              &   !gravitational constant
NA=6.02214076d-23,          &   !Avogadro Constant
kb=1.380658d-16,            &   !boltzmann constant
rsun=6.96d10,               &   !sun radius
rjupiter=7.14d9,            &   !jupiter radius
au=1.496d13,                &
km=1d5,                     &
year=3.15569d7,             &
day=86400d0,                &
t_hubble=4.55d17,           &   !hubble time
me=9.1093897d-28,           &   !electron mass
mp=1.6726231d-24,           &   !proton mass
c_light=2.9979d10,          &   !light speed
msun=1.99d33,               &   !sun mass
lsun=3.9d33,                &   !sun luminosity
mearth=5.976d27,            &   !earth mass
mjupiter=1.899d30,          &   !jupiter mass
h_planck=6.6260755d-27,     &   !planck constant
amu=1.660539040d-24,        &   !atomic mass unit
mh=1.6733d-24,              &   !hydrogen atom mass
mh2=3.3466d-24,             &   !hydrogen molecule mass
mhe=6.646481526d-24,        &   !helium mass
muh2=2.0158816d0,           &   !atomic weight of hydrogen molecule
muh=1.0076849d0,            &   !atomic weight of hydrogen
muhe=4.0026d0,              &   !atomic weight of helium
sigma_sb=5.67051d-5,        &   !Stefan-Boltzmann constant
ionh=2.18d-11,              &   !atomic hydrogen binding energy
dish=7.17d-12,              &   !molecular hydrogen dissociation energy
a_rad=7.5646d-15            !radiation energy density constant

real(8), parameter :: zero=0d0, fourth=0.25d0, half=0.5d0
real(8), parameter :: one=1d0, two=2d0, three=3d0, four=4d0
complex(8), parameter :: czero=(0.d0,0.d0)


interface initialize_global_parameters
    module procedure initialize_global_parameters_int,initialize_global_parameters_real,  &
    initialize_global_parameters_logical
end interface initialize_global_parameters

interface allocate_cell_data
    module procedure allocate_cell_data_1d,allocate_cell_data_nd
end interface allocate_cell_data

interface allocate_cell_data_block
    module procedure allocate_cell_data_block_1d,allocate_cell_data_block_nd
end interface allocate_cell_data_block

interface allocate_xsurface_data
    module procedure allocate_xsurface_data_1d,allocate_xsurface_data_nd
end interface allocate_xsurface_data

interface allocate_xsurface_data_block
    module procedure allocate_xsurface_data_block_1d,allocate_xsurface_data_block_nd
end interface allocate_xsurface_data_block

interface allocate_ysurface_data
    module procedure allocate_ysurface_data_1d,allocate_ysurface_data_nd
end interface allocate_ysurface_data

interface allocate_ysurface_data_block
    module procedure allocate_ysurface_data_block_1d,allocate_ysurface_data_block_nd
end interface allocate_ysurface_data_block

type axbsystem
    Vec x
    Vec b
    Mat A
    KSP ksp
end type axbsystem

type processor_variables
    logical :: amr_grow,smr_grow,processor_logical
    logical :: amr_derefine,smr_derefine
    logical :: petsc_decrease_rtol
    logical, dimension(:), allocatable :: mpi_logical,global_logical
    integer :: amr_ngrow,smr_ngrow
    integer :: amr_nderefine,smr_nderefine
    integer :: n_edge,n_corner
    integer, dimension(:), allocatable :: mpi_processor_oper
    integer, dimension(:), allocatable :: processor_operation_sum
    integer, dimension(:,:), allocatable :: grow_keys
    !each processor stores the sequence of the keys to a this list
    integer, dimension(:,:), allocatable :: derefine_keys
    integer, dimension(:,:), allocatable :: mpi_world_keys
    !rank0 processor stack all the keys to this list, then bcast it to the other processors
    real(8) :: collective_processor_result
    real(8) :: collective_result
    real(8), dimension(:), allocatable :: collective_array
end type processor_variables

type timedef
    real(8) :: t                            !current time
    real(8) :: t_output_initial             !the time of frame 1
    real(8) :: t_final                      !end program when t=tfinal
    real(8) :: dt_hydro                     !hydrodynamical timescale
    real(8) :: dt_processor                 !timestep of this processor
    real(8) :: dt_source                    !minimum of dt_thermal,dt_gravity,dt_radiation, and dt_hydro
    real(8) :: dt_thermal                   !thermal timescale
    real(8) :: dt_gravity                   !gravitational timescale
    real(8) :: dt_radiation                 !radiative timescale
    real(8) :: dt_frame                     !time interval between each frames
    real(8) :: dt_record                    !time interval to write a record
    real(8) :: t_next                       !next time of the output frame
    integer :: ntimestep
end type timedef

type scaledef
    !the default scales are cgs unit
    real(8) :: lscale                       !length scale
    real(8) :: timescale                    !timescale
    real(8) :: rscale                       !density scale
    real(8) :: tempscale                    !temperature scale
    real(8) :: vscale                       !velocity scale
    real(8) :: mscale                       !mass scale
end type scaledef

integer, parameter :: irho=1,ivx=2,ivy=3,ivz=4,ipres=5,itemp=6,iegv=7
integer, parameter :: imass=1,imomx=2,imomy=3,imomz=4,iE=5
integer, parameter :: pyro80=1,pyro70=2,pyro60=3,pyro50=4
integer, parameter :: ku_hlle=3,kl_hlle=3,ktotal=10
integer, parameter :: ku_lapack=2,kl_lapack=2,ktotal_lapack=7
integer, parameter :: transmissive=1,reflective=2,extrapolate=3,specified=9
integer, protected :: ioformat=2
integer, protected :: igravity=0    !self gravity=1, uniform gravity=2
integer, protected :: igeometry=0   !cartesian=0, cylindrical=1, spherical=2
integer, protected :: icooling=0
integer, protected :: iradiation=0
integer, protected :: nd
integer, protected :: nx
integer, protected :: ny
integer, protected :: nz
integer, protected :: llevel_max            !maxiumum logical level
integer, protected :: nx2
integer, protected :: n_hydro_guard
integer, protected :: xlb                   !x coordinate, lower boundary
integer, protected :: xub                   !x coordinate, upper boundary
integer, protected :: ylb
integer, protected :: yub
integer, protected :: nrefine_region
integer, protected :: max_refine_level
integer, protected :: np                    !number of processors
integer, protected :: rank
integer, protected :: restart_iframe
integer, protected :: blk_primitive_cell_size
integer, protected :: blk_cell_size
integer, protected :: blk_interface_size
integer, protected :: blk_xedge_size_hydro,blk_yedge_size_hydro,blk_corner_size_hydro
integer, public :: nmg
integer, public :: nradcell
integer, protected :: nray_static=32
logical, protected :: lradiation_advect=.false.
logical, protected :: lpassive=.false.
logical, protected :: lstiff_source
logical, protected :: restart
logical, public :: truncate_dt_hydro=.false.
logical, public :: resolve_rad_dt=.false.
logical, public :: viscous=.false.
logical, protected :: source=.false.
real(8), protected :: n_domain(4)
real(8), protected :: dxyz(3)
real(8), protected :: maw
real(8), protected :: gamma_gas
real(8), protected :: CFL
real(8), protected :: temp_Erad_boundary=100d0
real(8), protected :: fld_diffusion_boundary_value=0d0      
real(8), protected :: pfloor
real(8), protected :: g_uniform
real(8), public :: const_opacity=1d-2
real(8), public :: const_Erad=1d0
character(len=32), protected :: refine_type
character(len=16), protected :: east,south,west,north,se,sw,nw,ne
character(len=255), protected :: path_root
integer, protected :: blk_size_nx,blk_size_ny,blk_size_nz
integer, protected :: nx_blks,ny_blks,nz_blks,blk_xlb,blk_xub,blk_ylb,blk_yub

type ray2ddef
    !2d ray parameters
    real(8) :: I0,I                                 !the incoming intensity and the intensitiy solution
    real(8) :: source(3),chi(3),ds(2)
    real(8) :: pos(3)                               !xyz coordinate of the ray
end type ray2ddef

type blk_nb_coor
    !coordinates of immediate neighbouring cells, used in fld
    real(8), dimension(:), allocatable :: xl,xu,yl,yu
end type blk_nb_coor

type blockdef
    integer :: blk_id,level,static_level
    integer :: key(3)                       !(i,j,logical_level)
    integer :: nx_domain,ny_domain          !domain size of the logical node
    integer :: loc_type                     !determine the spatial type of the block
    real(8) :: dxyz(3)                      !coordinate difference of a cell
    !the lower limits of all three boundaries and the numerical size in x and y coordinates
    real(8) :: pos(3),n_size(2)
    real(8) :: Erad_pre,Erad_next,Erad_pre_int,Erad_next_int,sigma_rosseland_pre,sigma_rosseland_next
    real(8), dimension(:), allocatable :: Erad_xl,Erad_xu,Erad_yl,Erad_yu
    real(8), dimension(:), allocatable :: sigma_rosseland_xl,sigma_rosseland_xu,sigma_rosseland_yl,sigma_rosseland_yu
    real(8) :: dx_xl,dx_xu,dx_yl,dx_yu
    real(8), dimension(:), allocatable :: Erad_pre_mg,Erad_next_mg,sigma_rosseland_pre_mg,sigma_rosseland_next_mg
    real(8), dimension(:), allocatable :: dy
    real(8) :: dt_hydro_blk
    !Mesh points, for XDMF output
    real(8), dimension(:), allocatable :: mesh_x,mesh_y,curv
    !Surface quantities, values are on the interfaces
    real(8), dimension(:,:,:,:), allocatable :: xflux,yflux,flux_r,flux_theta,flux_phi
    real(8), dimension(:,:,:,:), allocatable :: xflux_predict,yflux_predict
    real(8), dimension(:,:,:), allocatable :: mass_flux_x,v1,v2,flux_r_amphi,flux_theta_amphi,flux_phi_amphi
    real(8), dimension(:,:,:), allocatable :: kx,ky,tau
    real(8), dimension(:,:,:), allocatable :: Fradx,Frady
    real(8), dimension(:,:,:,:), allocatable :: kx_mg,ky_mg
    real(8), dimension(:,:,:,:), allocatable :: Fradx_mg,Frady_mg
    real(8), dimension(:,:,:), allocatable :: surf1,surf2,lever_r,lever_theta
    real(8), dimension(:,:,:), allocatable :: gpotential_surf
    real(8), dimension(:,:,:), allocatable :: vis_tensor_xx,vis_tensor_yy,vis_tensor_xy,vis_tensor_yx
    !Volume quantities, values are at the center of the cells
    real(8), dimension(:,:,:,:), allocatable :: w,u,source,w1,w2,u1,u2
    real(8), dimension(:,:,:,:), allocatable :: predict_w,predict_u
    real(8), dimension(:,:,:,:), allocatable :: w0,u0
    real(8), dimension(:,:,:,:), allocatable :: w_xl,w_xr,w_yl,w_yr
    real(8), dimension(:,:,:,:), allocatable :: predict_fx,predict_fy
    real(8), dimension(:,:,:,:), allocatable :: dudt_x,dudt_y
    real(8), dimension(:,:,:,:), allocatable :: vol_center
    real(8), dimension(:,:,:), allocatable :: vol,lever,cs
    real(8), dimension(:,:,:), allocatable :: gpotential
    real(8), dimension(:,:,:), allocatable :: E_total
    real(8), dimension(:,:,:), allocatable :: vphi,vphi2,am_phi,am_phi2,am_theta,am_theta2,ekphi,ek,ek2
    real(8), dimension(:,:,:), allocatable :: vis_heat,vis_tensor_rphi,vis_tensor_thetaphi,vis_tensor_phitheta,vis_tensor_phiphi,divv,surf_phi
    !W primitive, U conservative, S source form. x,y,z (or i,j,k) and physical quantities
    !W physical quantities: 1=rho,2=vx,3=vy,4=vz,5=p
    !U physical quantities: 1=rho,2=rho vx,3=rho vy,4=rho vz,5=E
    real(8), dimension(:,:,:), allocatable :: egv,temp,cv,temp1,temp2,egv1,egv2,entropy
    real(8), dimension(:,:,:,:), allocatable :: xslp,yslp
    real(8), dimension(:,:,:), allocatable :: predict_temp,predict_egv
    real(8), dimension(:,:,:), allocatable :: temp0,egv0,temp_xl,temp_xr,temp_yl,temp_yr,egv_xl,egv_xr,egv_yl,egv_yr
    real(8), dimension(:,:,:), allocatable :: petsc_egv
    real(8), dimension(:,:,:), allocatable :: H2,HI,HII,electron
    real(8), dimension(:,:,:,:), allocatable :: s_geo1,s_geo2,s_grav
    real(8), dimension(:,:,:), allocatable :: Erad,Erad_int,feddington
    real(8), dimension(:,:,:), allocatable :: kappa_planck,kappa_rosseland                  !unit cm^2/g
    real(8), dimension(:,:,:), allocatable :: sigma_planck,sigma_rosseland
    real(8), dimension(:,:,:), allocatable :: kappa_sc,chi_sc,epsilon_sc
    real(8), dimension(:,:,:), allocatable :: rad_source
    real(8), dimension(:,:,:), allocatable :: gxl,gxr,gyl,gyr
    real(8), dimension(:,:,:,:), allocatable :: Erad_mg,feddington_mg
    real(8), dimension(:,:,:,:), allocatable :: a_mg,b_mg
    real(8), dimension(:,:,:,:), allocatable :: kappa_planck_mg,kappa_rosseland_mg
    real(8), dimension(:,:,:,:), allocatable :: sigma_planck_mg,sigma_rosseland_mg
    real(8), dimension(:), allocatable :: x_center,y_center,x_interface,y_interface
    real(8), dimension(:), allocatable :: stheta_center,ctheta_center,r_center,r_interface,theta_interface,theta_center
    real(8), dimension(:), allocatable :: stheta_interface,ctheta_interface
    logical :: derefine
    logical :: discontinuity
    logical :: on_processor
    logical :: reverse
    type(blk_nb_coor) :: nb_coor

    type(blockdef), pointer :: blk_pre,blk_next,blk_p,blk_xl,blk_xu,blk_xlyl,blk_xuyl,blk_xlyu,blk_xuyu
    type(blockneighbour), pointer :: nb_l,nb_r
    type(blockneighbour), pointer :: nb_n,nb_s,nb_e,nb_w,nb_ne,nb_nw,nb_se,nb_sw,nb_i
    type(ray2ddef), dimension(:,:,:), allocatable :: rays
end type blockdef

type blockneighbour
    integer :: cpu_rank         !for MPI
    type(blockdef), pointer :: blk
    type(blockneighbour), pointer :: blk_next       !two blocks are adjacent to one block
end type blockneighbour

type environmentdef
    !light data, parameters about the simulation
    real(8) :: h_ratio
    real(8) :: he_ratio
    real(8) :: nrho
    real(8) :: ntemp
    real(8) :: rhomin
    real(8) :: rhomax
    real(8) :: tmin
    real(8) :: tmax
    real(8), dimension(:), allocatable :: rho_eos,t_eos
    integer :: table_dim(2)
    integer :: meta_table_dim(3)
    logical :: lenforce_field_quantities
    logical :: firststep
end type environmentdef

type table1d
    real(8), dimension(:), allocatable :: t1d,xlist
    real(8) :: dx
    integer :: nd1
end type table1d

type table2d
    real(8), dimension(:,:), allocatable :: t2d          !values of the table
    real(8), dimension(:), allocatable :: xlist,ylist    !lists of sampled points
    real(8) :: dx,dy                                     !space of sampled points
    integer :: nd1,nd2                                  !sampled points
end type table2d

type uniform_table1d
    real(8) :: xmin,xmax,dx
end type uniform_table1d

type table_group
    type(table2d), dimension(:), allocatable :: table_eos
end type table_group

type particle
    real(8) :: mass
    real(8), dimension(3) :: xyz
end type particle

type stellar
    type(particle) :: core
    real(8) :: teff
    real(8) :: radius
    real(8) :: period
    real(8) :: spin
    real(8) :: luminosity               !blackbody part
    real(8) :: lum_atm                  !atmosphere luminosity
end type stellar

type(table_group), public :: tables
type(particle) :: single_star
type(stellar) :: central_star
logical, public :: EOS_test=.false.

real(8), public :: opacity_gas_rho_min=1d-18
real(8), public :: opacity_gas_rho_max=1d-5
real(8), public :: opacity_gas_pres_min
real(8), public :: opacity_gas_pres_max
real(8), public :: opacity_gas_t_min=200
real(8), public :: opacity_gas_t_max=40000
real(8), public :: rho_slope_refine_thresh,erad_slope_refine_thresh
real(8), dimension(:,:), allocatable, public :: temp_division
real(8), dimension(:,:), allocatable, public :: record_history
real(8), dimension(:,:), allocatable, public :: refine_region
real(8), dimension(:,:), allocatable, public :: rad_mg
integer, public :: rad_bound_type(4),hydro_bound_type(4)
integer :: record_i=0,record_chunk_i=0,record_array_size=0,record_chunk_size=100
integer, public :: metallicity=2                !1~-0.3 (metal poor), 2~0.0 (solar), 3~0.3 (metal rich)
integer, public :: nframe                       !number of frames to be generated
integer, public :: iframe                       !current frame number
integer, public :: record_length,j_record
integer, dimension(:), allocatable, public :: np_nblk
type(uniform_table1d) :: tab_rho_t
logical :: division_known
logical, public :: energy_conserve_formalism=.false.
logical, public :: rad_adv=.false.
real(8) :: timescale
real(8), dimension(:), allocatable, public :: dt_processors
real(8), dimension(:), allocatable, public :: t_target_output
real(8), dimension(:), allocatable, public :: theta_center,theta_interface
real(8), dimension(:), allocatable, public :: dr2dr3_ratio
real(8), public :: r_forbid
real(8), public :: petsc_rtol
real(8), public :: petsc_eos_rtol=1d-2
character(len=32), public :: out_dir='out'
character(len=128), public :: path_out
type(timedef), public :: time_sys
type(scaledef), public :: scale_sys
type(axbsystem), public :: axb_operator
!blocks related
integer, public :: nblk_total,nblk_processor
integer, public :: nblocks_tree
type(blockdef), public, pointer :: blk_head,blk_tail,blk_root,blk_processor_head,blk_processor_tail
type(blockdef), public, pointer :: blk_head_xl_bound,blk_head_xu_bound,blk_head_yl_bound,blk_head_yu_bound
type(environmentdef), public :: environment
type(processor_variables), public :: processor


procedure(condition_scalar), pointer :: boundcond_scalar,initialcond_scalar

contains

!prototype subroutines

subroutine collective_oper_sub1(blk,blk_result)
    !the prototype of collective operation, sub1
    !carry out operation on a block of a processor, save the result to blk_result
    type(blockdef), pointer :: blk
    real(8) :: blk_result
end subroutine collective_oper_sub1

subroutine collective_oper_sub2(blk_result_array)
    !the prototype of collective operation, sub2
    !each processor will summarize its result and save it to processor%collective_processor_result
    real(8), dimension(:), allocatable :: blk_result_array
end subroutine collective_oper_sub2

subroutine collective_oper_sub3()
    !the prototype of collective operation, sub3
    !rank0 will carry out this subroutine and save the result to processor%collective_result
    !then bcast the result to all other processors
end subroutine collective_oper_sub3

subroutine condition_hydro(pos,t,w,u,temp,egv,mark)
    !this prototype is used both to specify hydro conditions inside the computational domain
    !and at the boundaries, mark is used to distinguish different boundary locations or
    !some special cases
    real(8) :: pos(3),t,w(5),u(5),temp,egv
    integer, optional :: mark
end subroutine condition_hydro

subroutine condition_scalar(pos,t,scalar)
    !the prototype of subroutine that specifies a scalar
    real(8) :: pos(3),t,scalar
end subroutine condition_scalar

subroutine condition_rad(pos,t,Erad,Frad,feddington,mark)
    !the prototype of subroutine that specify Erad, Frad, and feddington at pos (xyz coordinate) and t
    real(8) :: pos(3),t,Frad(3),feddington
    real(8), dimension(:), allocatable :: Erad
    integer, optional :: mark
end subroutine condition_rad

subroutine condition_domain(pos,t,w,u,temp,egv,assign_value)
    !the prototype of subroutine that specify w, u, temp and egv at pos (xyz coordinate) and t
    real(8) :: pos(3),t,w(5),u(5),temp,egv
    logical :: assign_value
end subroutine condition_domain

subroutine condition_domainrad(pos,t,Erad,Frad,feddington,assign_value)
    !the prototype of subroutine that specify Erad, Frad and feddington at pos (xyz coordinate) and t
    real(8) :: pos(3),t,Erad,Frad(3),feddington
    logical :: assign_value
end subroutine condition_domainrad

subroutine extractarray(blk,ap)
    !the prototype of subroutine that extract an array from a block
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
end subroutine extractarray

subroutine blk_operator(blk)
    !the prototype of subroutine of blk operator, works for blk_traversal subroutine
    type(blockdef), pointer :: blk
end subroutine blk_operator

subroutine initialize_simulation_parameters()
    real(8) :: tfinal,v_boundary(4),h_ratio,he_ratio,rhomin,rhomax,tmin,tmax,     &
        refine_xmin,refine_xmax,refine_ymin,refine_ymax,r,r1,r2
    integer :: i,nrho,ntemp,level,ierr,nc,errcode
    namelist /meshinfo/ n_domain,tfinal,timescale,CFL,v_boundary,hydro_bound_type,rad_bound_type,    &
        nframe,refine_type,nrefine_region,max_refine_level,restart,restart_iframe
    namelist /phyinfo/ h_ratio,he_ratio
    namelist /global_parameters/ igravity,igeometry,icooling,iradiation,nd,nx,ny,    &
        blk_size_nx,blk_size_ny,maw,gamma_gas
    namelist /refinement/ refine_xmin,refine_xmax,refine_ymin,refine_ymax,level
    call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
    call mpi_comm_size(MPI_COMM_WORLD,np,ierr)
    call getcwd(path_root)
    east='east';south='south';west='west';north='north';se='se';sw='sw';nw='nw';ne='ne'
    allocate(np_nblk(np))
    if (rank==0) then
        open(unit=11,file='modules/problem/global.data',status='old',action='read')
        read(unit=11,nml=meshinfo)
        read(unit=11,nml=phyinfo)
        read(unit=11,nml=global_parameters)
        if (igeometry==1.or.igeometry==2) then
            n_domain(3:4)=n_domain(3:4)*pi
            print *,n_domain
        end if
    end if
    nz=1;blk_size_nz=1
    !the order should be consistent with the order of global.data
    call mpi_bcast(n_domain,4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(tfinal,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(timescale,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(CFL,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(v_boundary,4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hydro_bound_type,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rad_bound_type,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nframe,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(refine_type,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nrefine_region,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(max_refine_level,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(restart_iframe,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(h_ratio,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(he_ratio,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !call mpi_bcast(ioformat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(igravity,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(igeometry,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(icooling,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(iradiation,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(blk_size_nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(blk_size_ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(maw,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(gamma_gas,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    if (restart) then
        iframe=restart_iframe
    else
        iframe=0
    end if
    rho_slope_refine_thresh=0.3
    erad_slope_refine_thresh=0.02
    time_sys%ntimestep=0
    time_sys%t_final=tfinal*timescale
    time_sys%dt_record=0d0
    time_sys%dt_hydro=0d0
    if (refine_type=='none') then                                       !no mesh refinement
    else if (refine_type=='static'.or.refine_type=='mixed') then        !static mesh refinement
        allocate(refine_region(nrefine_region,5))
        do i=1,nrefine_region
            if (rank==0) then
                read(11,nml=refinement)
                refine_region(i,1)=refine_xmin
                refine_region(i,2)=refine_xmax
                refine_region(i,3)=refine_ymin
                refine_region(i,4)=refine_ymax
                refine_region(i,5)=level
            end if
        end do
        if (np>1) then
            nc=nrefine_region*5
            call mpi_bcast(refine_region,nc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        end if
    end if
    close(11)
    dxyz=(/(n_domain(2)-n_domain(1))/nx,(n_domain(4)-n_domain(3))/ny,0d0/)
    nrho=1201;ntemp=1001;rhomin=1d-20;rhomax=1d0;tmin=5d2;tmax=2d5
    environment%h_ratio=h_ratio
    environment%he_ratio=he_ratio
    environment%nrho=nrho
    environment%ntemp=ntemp
    environment%rhomin=rhomin
    environment%rhomax=rhomax
    environment%tmin=tmin
    environment%tmax=tmax
    environment%firststep=.true.
    if (igeometry==2) then
        allocate(theta_center(-1:ny/2),theta_interface(-2:ny/2))
        allocate(dr2dr3_ratio(-1:nx+2))
        do i=-2,ny/2
            theta_interface(i)=n_domain(3)+i*dxyz(2)
        end do
        do i=-1,ny/2
            theta_center(i)=n_domain(3)+(i-half)*dxyz(2)
        end do
        do i=-1,nx+2
            r1=n_domain(1)+(i-1)*dxyz(1)
            r2=r1+dxyz(1)
            dr2dr3_ratio(i)=(r2**3d0-r1**3d0)/(r2**2d0-r1**2d0)
        end do
    end if
    if (igravity/=0) then
        source=.true.
    end if
    allocate(dt_processors(np),processor%processor_operation_sum(np))
    allocate(processor%collective_array(np),processor%global_logical(np))
    call check_global_parameter_consistency()
    call initialize_mesh()
    call blk_size_parameters()
end subroutine initialize_simulation_parameters

subroutine blk_size_parameters()
    if (nd==1) then
        blk_primitive_cell_size=5*(blk_size_nx+2*n_hydro_guard)
        blk_cell_size=blk_size_nx+2*n_hydro_guard
        blk_interface_size=blk_cell_size-1
    else if (nd==2) then
        blk_xedge_size_hydro=2*blk_size_nx*5
        blk_yedge_size_hydro=2*blk_size_ny*5
        blk_corner_size_hydro=2*2*5
        blk_primitive_cell_size=5*(blk_size_nx+2*n_hydro_guard)*(blk_size_ny+2*n_hydro_guard)
        blk_cell_size=(blk_size_nx+2*n_hydro_guard)*(blk_size_ny+2*n_hydro_guard)
        blk_interface_size=(blk_size_nx+2*n_hydro_guard-1)*(blk_size_ny+2*n_hydro_guard)
    end if
end subroutine blk_size_parameters

subroutine create_table_1d(m,table)
    integer :: m
    type(table1d) :: table
    allocate(table%t1d(m),table%xlist(m))
    table%nd1=m
end subroutine create_table_1d

subroutine create_table_2d(m,n,table)
    integer :: m,n
    type(table2d) :: table
    allocate(table%t2d(m,n),table%xlist(m),table%ylist(n))
    table%nd1=m
    table%nd2=n
end subroutine create_table_2d

subroutine calculate_llevel_max()
    integer :: nxblock,nyblock
    character(len=128) :: alert
    if (nd==1) then
        if (mod(nx,blk_size_nx)/=0) then
            alert='block size is not a divider of domain size'
            call abort_achilles(alert)
        end if
        nxblock=nx/blk_size_nx
        llevel_max=0
        do while (nxblock/=1)
            nxblock=ceiling(real(nxblock)/2)
            llevel_max=llevel_max+1
        end do
    else if (nd==2) then
        if (mod(nx,blk_size_nx)/=0.or.mod(ny,blk_size_ny)/=0) then
            alert='block size is not a divider of domain size'
            call abort_achilles(alert)
        end if
        nxblock=nx/blk_size_nx
        nyblock=ny/blk_size_ny
        llevel_max=0
        do while (nxblock/=1.or.nyblock/=1)
            nxblock=ceiling(real(nxblock)/2)
            nyblock=ceiling(real(nyblock)/2)
            llevel_max=llevel_max+1
        end do
    end if
end subroutine calculate_llevel_max

subroutine block_key_to_pointer(key,blk)
    !given a block key, return the pointer to the key (blk)
    !if such a block does not exist, return the *closest* parent key
    type(blockdef), pointer :: blk
    integer :: key(3),i,llevel,idx,idy
    integer, allocatable :: seq(:)
    allocate(seq(key(3)));seq=-1
    if (nd==1) then
        i=1;llevel=key(3);idx=key(1)
        do while (llevel>0)
            seq(i)=mod(idx,2)
            i=i+1
            llevel=llevel-1
            idx=ceiling(real(idx)/2)
        end do
        blk=>blk_root
        do i=key(3),1,-1
            if (seq(i)==1) then
                if (associated(blk%blk_xl)) then
                    blk=>blk%blk_xl
                else
                    !no deeper levels
                    exit
                end if
            else if (seq(i)==0) then
                if (associated(blk%blk_xu)) then
                    blk=>blk%blk_xu
                else
                    !no deeper levels
                    exit
                end if
            end if
        end do
    else if (nd==2) then
        i=1;llevel=key(3);idx=key(1);idy=key(2)
        do while (llevel>0)
            if (mod(idx,2)==1.and.mod(idy,2)==1) then
                seq(i)=1
            else if (mod(idx,2)==0.and.mod(idy,2)==1) then
                seq(i)=2
            else if (mod(idx,2)==1.and.mod(idy,2)==0) then
                seq(i)=3
            else if (mod(idx,2)==0.and.mod(idy,2)==0) then
                seq(i)=4
            else
                print *, 'wrong key_to_pointer'
            end if
            i=i+1
            llevel=llevel-1
            idx=ceiling(real(idx)/2)
            idy=ceiling(real(idy)/2)
        end do
        blk=>blk_root
        do i=key(3),1,-1
            select case (seq(i))
            case (1)
                if (associated(blk%blk_xlyl)) then
                blk=>blk%blk_xlyl
                else
                exit
                end if
            case (2)
                if (associated(blk%blk_xuyl)) then
                blk=>blk%blk_xuyl
                else
                exit
                end if
            case (3)
                if (associated(blk%blk_xlyu)) then
                blk=>blk%blk_xlyu
                else
                exit
                end if
            case (4)
                if (associated(blk%blk_xuyu)) then
                blk=>blk%blk_xuyu
                else
                exit
                end if
                case default
                print *, 'wrong key_to_pointer seq'
            end select
        end do
    end if
    deallocate(seq)
end subroutine block_key_to_pointer

function key_xl_neighbour(key)
    !return the key of the neighbour block on x lower side
    integer :: key(3),key_xl_neighbour(3)
    if (nd==1) then
        key_xl_neighbour(1)=key(1)-1
        key_xl_neighbour(2:3)=key(2:3)
    else if (nd==2) then
    end if
end function key_xl_neighbour

function key_xu_neighbour(key)
    !return the key of the neighbour block on x upper side
    integer :: key(3),key_xu_neighbour(3)
    if (nd==1) then
        key_xu_neighbour(1)=key(1)+1
        key_xu_neighbour(2:3)=key(2:3)
    else if (nd==2) then
    end if
end function key_xu_neighbour

function key_parent(key)
    integer :: key(3),key_parent(3)
    if (nd==1) then
        key_parent(1)=ceiling(real(key(1))/2)
        key_parent(2)=key(2)
        key_parent(3)=key(3)-1
    else if (nd==2) then
    end if
end function key_parent

    subroutine logical_level_order_to_block_key(level_order,key)
integer :: key(3),level_order(llevel_max)
    end subroutine logical_level_order_to_block_key

subroutine build_logical_tree()
    !suppose we have square blocks, i.e., blk_size_nx=blk_size_ny
    !decompose the whole domain into logical domain until the smallest logical domain
    !equals the coarsest block
    integer :: i
    if (nd==1) then
        i=0
        call new_blank_block(blk_root)
        blk_root%key=(/1,1,0/)
        blk_root%nx_domain=nx
        call new_binary_block_recursive(blk_root,i)
    else if (nd==2) then
        i=0
        call new_blank_block(blk_root)
        blk_root%key=(/1,1,0/)
        blk_root%nx_domain=nx
        blk_root%ny_domain=ny
        call new_quad_block_recursive(blk_root,i)
    end if
end subroutine build_logical_tree

subroutine new_binary_block_recursive(blk,i)
    type(blockdef), pointer :: blk
    integer :: i,j,level_size_max,level_size,level_size_next
    character(len=128) :: alert
    level_size_max=blk_size_nx*2**(llevel_max-1)
    level_size=blk_size_nx*2**(llevel_max-i)
    j=i
    if (j<=llevel_max.and.level_size/=blk_size_nx) then
        level_size_next=level_size/2
        j=j+1
        call new_blank_block(blk%blk_xl)
        blk%blk_xl%nx_domain=level_size_next
        blk%blk_xl%key=(/2*blk%key(1)-1,1,j/)
        blk%blk_xl%blk_p=>blk
        call new_binary_block_recursive(blk%blk_xl,j)
        if (blk%nx_domain>level_size_next) then
            call new_blank_block(blk%blk_xu)
            blk%blk_xu%nx_domain=min(blk%nx_domain-level_size_next,level_size_next)
            blk%blk_xu%key=(/2*blk%key(1),1,j/)
            blk%blk_xu%blk_p=>blk
            blk%blk_xl%blk_next=>blk%blk_xu
            call new_binary_block_recursive(blk%blk_xu,j)
        end if
    end if
end subroutine new_binary_block_recursive

subroutine new_quad_block_recursive(blk,i)
    !2d Z-ordering xlyl->xuyl->xlyu->xuyu
    type(blockdef), pointer :: blk
    integer :: i,j,level_size_max,level_size,level_size_next
    character(len=128) :: alert
    level_size_max=blk_size_nx*2**(llevel_max-1)
    level_size=blk_size_nx*2**(llevel_max-i)
    j=i
    if (j<=llevel_max.and.level_size/=blk_size_nx) then
        level_size_next=level_size/2
        j=j+1
        call new_blank_block(blk%blk_xlyl)
        blk%blk_xlyl%nx_domain=min(level_size_next,blk%nx_domain)
        blk%blk_xlyl%ny_domain=min(level_size_next,blk%ny_domain)
        blk%blk_xlyl%key=(/2*blk%key(1)-1,2*blk%key(2)-1,j/)
        blk%blk_xlyl%blk_p=>blk
        call new_quad_block_recursive(blk%blk_xlyl,j)
        if (blk%nx_domain>level_size_next) then
            if (blk%ny_domain>level_size_next) then
                call new_blank_block(blk%blk_xuyl)
                call new_blank_block(blk%blk_xlyu)
                call new_blank_block(blk%blk_xuyu)
                blk%blk_xuyl%nx_domain=min(blk%nx_domain-level_size_next,level_size_next)
                blk%blk_xuyl%ny_domain=level_size_next
                blk%blk_xuyl%key=(/2*blk%key(1),2*blk%key(2)-1,j/)
                blk%blk_xuyl%blk_p=>blk
                blk%blk_xlyl%blk_next=>blk%blk_xuyl
                blk%blk_xlyu%nx_domain=level_size_next
                blk%blk_xlyu%ny_domain=min(blk%ny_domain-level_size_next,level_size_next)
                blk%blk_xlyu%key=(/2*blk%key(1)-1,2*blk%key(2),j/)
                blk%blk_xlyu%blk_p=>blk
                blk%blk_xuyl%blk_next=>blk%blk_xlyu
                blk%blk_xuyu%nx_domain=min(blk%nx_domain-level_size_next,level_size_next)
                blk%blk_xuyu%ny_domain=min(blk%ny_domain-level_size_next,level_size_next)
                blk%blk_xuyu%key=(/2*blk%key(1),2*blk%key(2),j/)
                blk%blk_xuyu%blk_p=>blk
                blk%blk_xlyu%blk_next=>blk%blk_xuyu
                call new_quad_block_recursive(blk%blk_xuyl,j)
                call new_quad_block_recursive(blk%blk_xlyu,j)
                call new_quad_block_recursive(blk%blk_xuyu,j)
            else
                call new_blank_block(blk%blk_xuyl)
                blk%blk_xuyl%nx_domain=min(blk%nx_domain-level_size_next,level_size_next)
                blk%blk_xuyl%ny_domain=blk%ny_domain
                blk%blk_xuyl%key=(/2*blk%key(1),2*blk%key(2)-1,j/)
                blk%blk_xuyl%blk_p=>blk
                blk%blk_xlyl%blk_next=>blk%blk_xuyl
                call new_quad_block_recursive(blk%blk_xuyl,j)
            end if
        else
            if (blk%ny_domain>level_size_next) then
                call new_blank_block(blk%blk_xlyu)
                blk%blk_xlyu%nx_domain=blk%nx_domain
                blk%blk_xlyu%ny_domain=min(blk%ny_domain-level_size_next,level_size_next)
                blk%blk_xlyu%key=(/2*blk%key(1)-1,2*blk%key(2),j/)
                blk%blk_xlyu%blk_p=>blk
                blk%blk_xlyl%blk_next=>blk%blk_xlyu
                call new_quad_block_recursive(blk%blk_xlyu,j)
            else
                !there is only one logical sub-domain in this level of domain
            end if
        end if
    end if
end subroutine new_quad_block_recursive

subroutine display_logical_tree()
    integer :: i
    if (nd==1) then
        i=0
        call display_logical_tree_node(blk_root,i)
    else if (nd==2) then
        i=0
        call display_logical_tree_node(blk_root,i)
    end if
end subroutine display_logical_tree

subroutine display_logical_tree_node(blk,i)
    type(blockdef), pointer :: blk
    integer :: i,j
    j=i
    if (nd==1) then
        if (associated(blk)) then               !this level exist
            print *,blk%key,associated(blk%blk_xl),associated(blk%blk_xu),associated(blk%blk_p)
            j=j+1
            if (associated(blk%blk_xl)) then
                call display_logical_tree_node(blk%blk_xl,j)
            end if
            if (associated(blk%blk_xu)) then
                call display_logical_tree_node(blk%blk_xu,j)
            end if
        end if
    else if (nd==2) then
        if (associated(blk)) then               !this level exist
            print *,blk%key,associated(blk%blk_xlyl),associated(blk%blk_xuyl),associated(blk%blk_xlyu),associated(blk%blk_xuyu)
            j=j+1
            if (associated(blk%blk_xlyl)) then
                call display_logical_tree_node(blk%blk_xlyl,j)
            end if
            if (associated(blk%blk_xuyl)) then
                call display_logical_tree_node(blk%blk_xuyl,j)
            end if
            if (associated(blk%blk_xlyu)) then
                call display_logical_tree_node(blk%blk_xlyu,j)
            end if
            if (associated(blk%blk_xuyu)) then
                call display_logical_tree_node(blk%blk_xuyu,j)
            end if
        end if
    end if
end subroutine display_logical_tree_node

subroutine find_the_first_node(blk1,blk2)
    !find the first node of branch blk1, store in blk2
    type(blockdef), pointer :: blk1,blk2
    if (nd==1) then
    if (associated(blk1%blk_xl)) then
call find_the_first_node(blk1%blk_xl,blk2)
    else
    blk2=>blk1
    end if
    else if (nd==2) then
    if (associated(blk1%blk_xlyl)) then
call find_the_first_node(blk1%blk_xlyl,blk2)
    else
    blk2=>blk1
    end if
    end if
    end subroutine find_the_first_node

subroutine find_the_successor_node(blk1,blk2)
    type(blockdef), pointer :: blk1,blk2
    if (nd==1) then
        if (associated(blk1%blk_next)) then             !the next logical block is in this level
            blk2=>blk1%blk_next
        else
            blk2=>blk1%blk_p                            !go to the parent levels to find the next logical block
            do while (.not.associated(blk2%blk_next))
                if (associated(blk2,blk_root)) then
                blk2=>blk_root
                goto 10
                end if
                blk2=>blk2%blk_p
            end do
            call find_the_first_node(blk2%blk_next,blk2)
10      end if
    else if (nd==2) then
        if (associated(blk1%blk_next)) then             !the next logical block is in this level
            blk2=>blk1%blk_next
        else                                            !go to the parent levels to find the next logical block
            blk2=>blk1%blk_p
        do while (.not.associated(blk2%blk_next))
            if (associated(blk2,blk_root)) then     !the end of the tree
                blk2=>blk_root
                goto 20
            end if
            blk2=>blk2%blk_p
        end do
        call find_the_first_node(blk2%blk_next,blk2)
20      end if
    end if
end subroutine find_the_successor_node

subroutine calculate_domain_block_parameters()
    nx_blks=nx/blk_size_nx
    ny_blks=ny/blk_size_ny
    nz_blks=nz/blk_size_nz
    blk_xlb=1-n_hydro_guard
    blk_xub=blk_size_nx+n_hydro_guard
end subroutine calculate_domain_block_parameters

subroutine new_blank_block(blk)
    type(blockdef), pointer :: blk
    allocate(blk)
    blk%blk_id=-1
    blk%level=-1
    blk%static_level=0
    nullify(blk%blk_pre,blk%blk_next,blk%blk_p,blk%blk_xl,blk%blk_xu,   &
        blk%blk_xlyl,blk%blk_xuyl,blk%blk_xlyu,blk%blk_xuyu)
    nullify(blk%nb_n,blk%nb_s,blk%nb_e,blk%nb_w,blk%nb_ne,blk%nb_nw,blk%nb_se,blk%nb_sw,blk%nb_i)
    nullify(blk%nb_l,blk%nb_r)
    blk%derefine=.false.
    blk%discontinuity=.false.
    blk%on_processor=.false.
end subroutine new_blank_block

subroutine new_block(blk,level,coords,n_size)
    type(blockdef), pointer :: blk
    integer :: id,level,i,j,k
    real(8) :: coords(2),n_size(2)
    blk%dxyz=dxyz/2d0**level
    blk%level=level
    blk%pos(1:2)=coords
    blk%pos(3)=0
    blk%n_size=n_size
    if (nd==1) then
        allocate(blk%mesh_x(blk_size_nx+1))
        do i=1,blk_size_nx+1
            blk%mesh_x(i)=coords(1)+(i-1)*blk%dxyz(1)
        end do
    else if (nd==2) then
        allocate(blk%mesh_x(blk_size_nx+1),blk%mesh_y(blk_size_ny+1))
        do i=1,blk_size_nx+1
            blk%mesh_x(i)=coords(1)+(i-1)*blk%dxyz(1)
        end do
        do i=1,blk_size_ny+1
            blk%mesh_y(i)=coords(2)+(i-1)*blk%dxyz(2)
        end do
    end if
end subroutine new_block

subroutine build_linked_base_blocks()
    !currently, all processors have the same linked list.
    !the linked list has basic geometric information
    type(blockdef), pointer :: blk_builder1,blk_builder2
    integer :: i,j,k,id,level
    real(8) :: coords(2),n_size(2)
    level=0
    call calculate_domain_block_parameters()
    if (nd==1) then
        call find_the_first_node(blk_root,blk_head)
        coords=(/n_domain(1),0d0/)
        n_size=(/(n_domain(2)-n_domain(1))/nx_blks,0d0/)
        call new_block(blk_head,level,coords,n_size)
        blk_builder1=>blk_head
        do while (associated(blk_builder1))
            call find_the_successor_node(blk_builder1,blk_builder2)
            if (.not.associated(blk_builder2,blk_root)) then
                blk_builder1%blk_next=>blk_builder2
                blk_builder2%blk_pre=>blk_builder1
                blk_builder1=>blk_builder2
                coords(1)=n_domain(1)+(n_domain(2)-n_domain(1))/nx_blks*(blk_builder1%key(1)-1)
                coords(2)=0d0
                n_size(1)=(n_domain(2)-n_domain(1))/nx_blks
                n_size(2)=0d0
                call new_block(blk_builder1,level,coords,n_size)
            else
                exit
            end if
        end do
        blk_tail=>blk_builder1
        nullify(blk_builder1,blk_builder2)
    else if (nd==2) then
        blk_ylb=1-n_hydro_guard
        blk_yub=blk_size_ny+n_hydro_guard
        coords(1)=n_domain(1)
        coords(2)=n_domain(3)
        call find_the_first_node(blk_root,blk_head)
        coords=(/n_domain(1),n_domain(3)/)
        n_size=(/(n_domain(2)-n_domain(1))/nx_blks,(n_domain(4)-n_domain(3))/ny_blks/)
        call new_block(blk_head,level,coords,n_size)
        blk_builder1=>blk_head
        do while (associated(blk_builder1))
            call find_the_successor_node(blk_builder1,blk_builder2)
            if (.not.associated(blk_builder2,blk_root)) then
                blk_builder1%blk_next=>blk_builder2
                blk_builder2%blk_pre=>blk_builder1
                blk_builder1=>blk_builder2
                coords(1)=n_domain(1)+(n_domain(2)-n_domain(1))/nx_blks*(blk_builder1%key(1)-1)
                coords(2)=n_domain(3)+(n_domain(4)-n_domain(3))/ny_blks*(blk_builder1%key(2)-1)
                n_size(1)=(n_domain(2)-n_domain(1))/nx_blks
                n_size(2)=(n_domain(4)-n_domain(3))/ny_blks
                call new_block(blk_builder1,level,coords,n_size)
            else
                exit
            end if
        end do
        blk_tail=>blk_builder1
        nblk_total=blk_tail%blk_id
        nullify(blk_builder1,blk_builder2)
    end if
end subroutine build_linked_base_blocks

subroutine new_neighbour_block(blk_nb)
    type(blockneighbour), pointer :: blk_nb
    allocate(blk_nb)
    nullify(blk_nb%blk_next)
end subroutine new_neighbour_block

subroutine display_blocks()
    type(blockdef), pointer :: blk
    nullify(blk)
    blk=>blk_processor_head
    do while (associated(blk))
        print *,blk%blk_id,blk%blk_pre%blk_id,blk%blk_next%blk_id,blk%key,blk%level,blk%static_level,blk%loc_type
        !write (*,'(8ES15.5E2)')blk%mesh_x
        blk=>blk%blk_next
    end do
end subroutine display_blocks

subroutine display_block_data(blk)
    !including the buffer zones
    type(blockdef), pointer :: blk
    character(len=32) :: s1,fmt1
    integer :: nbx,nby,i
    if (nd==1) then
    else if (nd==2) then
        nbx=blk_xub-blk_xlb+1
        nby=blk_yub-blk_ylb+1
        write(s1,fmt='(I5)') nby
        fmt1=trim("(")//trim(adjustl(s1))//trim("ES11.2E2)")
        print *,'density'
        do i=blk_xlb,blk_xub
            write(unit=6,fmt=fmt1) blk%w(i,blk_ylb:blk_yub,1,1)
        end do
        print *,'vx'
        do i=blk_xlb,blk_xub
            write(unit=6,fmt=fmt1) blk%w(i,blk_ylb:blk_yub,1,2)
        end do
        print *,'vy'
        do i=blk_xlb,blk_xub
            write(unit=6,fmt=fmt1) blk%w(i,blk_ylb:blk_yub,1,3)
        end do
        print *,'pressure'
        do i=blk_xlb,blk_xub
            write(unit=6,fmt=fmt1) blk%w(i,blk_ylb:blk_yub,1,5)
        end do
        print *,'egv'
        do i=blk_xlb,blk_xub
            write(unit=6,fmt=fmt1) blk%egv(i,blk_ylb:blk_yub,1)
        end do
    end if
end subroutine display_block_data

subroutine check_leaves()
    type(blockdef), pointer :: blk
    real(8) :: r
    nullify(blk)
    blk=>blk_head
    do while (associated(blk))
        if (dxyz(1)/blk%dxyz(1)/=2**blk%level) then
            print *,'level wrong',blk%key,blk%level,dxyz(1)/blk%dxyz(1)
            stop
        end if
        if (associated(blk%blk_pre)) then
            if (abs(blk%blk_pre%level-blk%level)>1) then
                print *,'1,check linked blocks, pre',blk%key,blk%blk_pre%key
                stop
            end if
            r=blk%blk_pre%dxyz(1)/blk%dxyz(1)
            if (r>2.or.r<0.5) then
                print *,'2,check linked blocks, pre',blk%key,blk%blk_pre%key,r
                stop
            end if
        end if
        if (associated(blk%blk_next)) then
            if (abs(blk%blk_next%level-blk%level)>1) then
                print *,'1,check linked blocks, next',blk%key,blk%blk_next%key
                stop
            end if
            r=blk%blk_next%dxyz(1)/blk%dxyz(1)
            if (r>2.or.r<0.5) then
                print *,'2,check linked blocks, next',blk%key,blk%blk_next%key,r
                stop
            end if
        end if
        blk=>blk%blk_next
    end do
end subroutine check_leaves

subroutine print_block_quantity(s)
    type(blockdef), pointer :: blk
    character(len=*) :: s
    blk=>blk_head
    print *,trim(s),time_sys%t,time_sys%ntimestep
    do while (associated(blk))
        select case (trim(s))
        case ('rho')
            write(*,'(2I8,8ES15.4E2)') blk%blk_id,blk%level,blk%w(1:blk_size_nx,1,1,1)
        case ('vx')
            write(*,'(2I8,8ES15.4E2)') blk%blk_id,blk%level,blk%w(1:blk_size_nx,1,1,2)
        case ('p')
            write(*,'(2I8,8ES15.4E2)') blk%blk_id,blk%level,blk%w(1:blk_size_nx,1,1,5)
        case ('temp')
            write(*,'(2I8,8ES15.4E2)') blk%blk_id,blk%level,blk%temp(1:blk_size_nx,1,1)
        end select
        blk=>blk%blk_next
    end do
end subroutine print_block_quantity

subroutine refine_static_tree_and_blocks()
    type(blockdef), pointer :: blk
    real(8) :: zone(4)
    integer :: i,static_level
    blk=>blk_head
    if (refine_type=='static'.or.refine_type=='mixed') then
        if (nd==1) then
            do i=1,nrefine_region
                zone=refine_region(i,1:4)
                static_level=int(refine_region(i,5))
                do while (associated(blk))
                    if (block_zone_overlap(blk,zone)) then
                        call grow_static_tree_node(blk,zone,static_level)
                    end if
                    blk=>blk%blk_next
                end do
            end do
            do i=1,max_refine_level-1
                blk=>blk_head
                do while (associated(blk))
                    call add_static_buffer_blocks(blk)
                    blk=>blk%blk_next
                end do
            end do
        else if (nd==2) then
            !no smr first
        end if
    end if
    nullify(blk)
end subroutine refine_static_tree_and_blocks

subroutine grow_static_tree_node(blk,zone,static_level)
    !recursively create deeper levels of blocks on all processors
    !does not allocate heavy data
    type(blockdef), pointer :: blk
    real(8) :: coords(2),n_size(2),zone(4)
    integer :: id,level,static_level
    if (nd==1) then
        call grow_tree_node_1d(blk)
        level=blk%level+1
        blk%blk_xl%static_level=level
        blk%blk_xu%static_level=level
        if (level<static_level) then
            if (block_zone_overlap(blk%blk_xl,zone)) then
                call grow_static_tree_node(blk%blk_xl,zone,static_level)
            end if
            if (block_zone_overlap(blk%blk_xu,zone)) then
                call grow_static_tree_node(blk%blk_xu,zone,static_level)
            end if
        end if
    else if (nd==2) then
    end if
end subroutine grow_static_tree_node

subroutine grow_tree_node_1d(blk)
    !grow a binary tree node and link
    !maintain the linked list
    !no heavy data manipulation
    type(blockdef), pointer :: blk
    real(8) :: coords(2),n_size(2)
    integer :: level
    call new_blank_block(blk%blk_xl)
    call new_blank_block(blk%blk_xu)
    blk%blk_xl%key=(/2*blk%key(1)-1,1,blk%key(3)+1/)
    blk%blk_xu%key=(/2*blk%key(1),1,blk%key(3)+1/)
    if (associated(blk,blk_processor_head)) then
        blk_processor_head=>blk%blk_xl
    end if
    if (associated(blk,blk_processor_tail)) then
        blk_processor_tail=>blk%blk_xu
    end if
    if (associated(blk%blk_pre)) then
        blk%blk_pre%blk_next=>blk%blk_xl
        blk%blk_xl%blk_pre=>blk%blk_pre
    else            !this is the head block
        blk_head=>blk%blk_xl
    end if
    if (associated(blk%blk_next)) then
        blk%blk_next%blk_pre=>blk%blk_xu
        blk%blk_xu%blk_next=>blk%blk_next
    else            !this is the tail block
        blk_tail=>blk%blk_xu
    end if
    blk%blk_xu%blk_pre=>blk%blk_xl
    blk%blk_xl%blk_next=>blk%blk_xu
    blk%blk_xl%blk_p=>blk
    blk%blk_xu%blk_p=>blk
    level=blk%level+1
    coords=blk%pos(1:2)
    n_size=blk%n_size/2
    call new_block(blk%blk_xl,level,coords,n_size)
    coords=blk%pos(1:2)+blk%n_size/2
    call new_block(blk%blk_xu,level,coords,n_size)
end subroutine grow_tree_node_1d

subroutine grow_tree_node_1d_conditional(blk)
    !if blk already has child blocks, do nothing, else grow the tree
    type(blockdef), pointer :: blk
    if (associated(blk%blk_xl)) then
        !do nothing
    else
        call grow_tree_node_1d(blk)
        call allocate_block_heavy_data(blk%blk_xl)
        call allocate_block_heavy_data(blk%blk_xu)
    end if
end subroutine grow_tree_node_1d_conditional

subroutine trim_tree_node(blk)
    type(blockdef), pointer :: blk
end subroutine trim_tree_node

function block_zone_overlap(blk,zone)
    logical :: block_zone_overlap
    type(blockdef), pointer :: blk
    real(8) :: zone(4),xmin,xmax,ymin,ymax
    if (nd==1) then
        xmin=blk%pos(1)
        xmax=blk%pos(1)+blk%n_size(1)
        if (xmin<=zone(1)) then
            if (xmax<=zone(1)) then
                block_zone_overlap=.false.
            else
                block_zone_overlap=.true.
            end if
        else
            if (xmin<=zone(2)) then
                block_zone_overlap=.true.
            else
                block_zone_overlap=.false.
            end if
        end if
    else if (nd==2) then
    end if
end function block_zone_overlap

subroutine add_static_buffer_blocks(blk)
    type(blockdef), pointer :: blk,blk_temp
    integer :: level
    if (nd==1) then
        if (associated(blk)) then
            level=blk%level
            if (associated(blk%blk_pre)) then
                blk_temp=>blk%blk_pre
                do while (blk_temp%level+1<level)
                    call grow_tree_node_1d(blk_temp)
                    blk_temp=>blk%blk_pre
                end do
            end if
            if (associated(blk%blk_next)) then
                blk_temp=>blk%blk_next
                do while (blk_temp%level+1<level)
                    call grow_tree_node_1d(blk_temp)
                    blk_temp=>blk%blk_next
                end do
            end if
        end if
    else if (nd==2) then
    end if
end subroutine add_static_buffer_blocks

subroutine identify_block_spatial_type()
    !       location type
    !           1   2   3
    !           4   5   6
    !           7   8   9
    type(blockdef), pointer :: blk
    integer :: i,key_x,key_y
    if (nd==1) then
        blk=>blk_head
        do i=1,nblk_total
            if (associated(blk,blk_head)) then
                blk%loc_type=1
            else if (associated(blk,blk_tail)) then
                blk%loc_type=3
            else
                blk%loc_type=2
            end if
            if (i/=nblk_total) then
                blk=>blk%blk_next
            end if
        end do
        nullify(blk)
    else if (nd==2) then
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            key_x=blk%key(1)
            key_y=blk%key(2)
            if (key_x==1) then
                if (key_y==1) then
                    blk%loc_type=7
                else if (key_y==ny_blks) then
                    blk%loc_type=1
                else
                    blk%loc_type=4
                end if
            else if (key_x==nx_blks) then
                if (key_y==1) then
                    blk%loc_type=9
                else if (key_y==ny_blks) then
                    blk%loc_type=3
                else
                    blk%loc_type=6
                end if
            else
                if (key_y==1) then
                    blk%loc_type=8
                else if (key_y==ny_blks) then
                    blk%loc_type=2
                else
                    blk%loc_type=5
                end if
            end if
            if (i/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
        nullify(blk)
    else
    end if
end subroutine identify_block_spatial_type

function is_xl_block(blk)
    type(blockdef), pointer :: blk
    logical :: is_xl_block
    if (blk%key(1)==1) then
        is_xl_block=.true.
    else
        is_xl_block=.false.
    end if
end function is_xl_block

function is_xu_block(blk)
    type(blockdef), pointer :: blk
    logical :: is_xu_block
    if (blk%key(1)==nx_blks*2**(blk%key(3)-llevel_max)) then
        is_xu_block=.true.
    else
        is_xu_block=.false.
    end if
end function is_xu_block

subroutine allocate_block_heavy_data(blk)
    type(blockdef), pointer :: blk
    call initialize_passive_quantities_block(blk)
    call initialize_geometry_block(blk)
    call initialize_hydro_block(blk)
    call initialize_eos_block(blk)
    call initialize_cooling_block(blk)
    call initialize_radiation_block(blk)
    call initialize_source_block(blk)
    call initialize_viscosity_block(blk)
end subroutine allocate_block_heavy_data

subroutine allocate_all_block_heavy_data()
    type(blockdef), pointer :: blk
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        call allocate_block_heavy_data(blk)
        blk=>blk%blk_next
    end do
    nullify(blk)
end subroutine allocate_all_block_heavy_data

subroutine renumber_domain_blocks()
    type(blockdef), pointer :: blk
    integer :: id
    id=1
    blk=>blk_head
    do while (associated(blk))
        blk_tail=>blk
        blk%blk_id=id
        id=id+1
        blk=>blk%blk_next
    end do
    nblk_total=blk_tail%blk_id
    nullify(blk)
end subroutine renumber_domain_blocks

subroutine destroy_block(blk)
    type(blockdef), pointer :: blk
    deallocate(blk)
end subroutine destroy_block

subroutine release_block(blk)
    !release the block heavy data, however, we actually do a substitution
    !the substituting block will inherit the structure, location, linked list and
    !neighbours of the old block. But not the heavy data.
    type(blockdef), pointer :: blk,blk_temp
    call destroy_block(blk)
end subroutine release_block

subroutine initialize_passive_quantities_block(blk)
    type(blockdef), pointer :: blk
    if (iradiation/=0) then
        call allocate_xsurface_data_block(blk%mass_flux_x)
    end if
end subroutine initialize_passive_quantities_block

subroutine initialize_geometry_block(blk)
    !surf1 is the surface area (line length in 2d) that is perpendicular to the 1st coordinate, same for surf2
    type(blockdef), pointer :: blk
    real(8) :: coords(3),dx,dxyz(3),r,theta,theta1,theta2,r1,r2,v1,v2,v,pos(3),r_val,sin_val,cos_val,theta_remap
    integer :: i,j,ijk(3),k
    if (nd==1) then
        coords=blk%pos
        dx=blk%dxyz(1)
        call allocate_cell_data_block(blk%vol)
        allocate(blk%x_center(blk_xlb:blk_xub),blk%x_interface(blk_xlb-1:blk_xub),blk%surf1(blk_xlb-1:blk_xub,1,1))
        if (igeometry==0) then
            do i=blk_xlb,blk_xub
                blk%x_center(i)=coords(1)+(i-half)*dx
                blk%vol(i,1,1)=dx
            end do
            do i=blk_xlb-1,blk_xub
                blk%surf1=1d0
                blk%x_interface(i)=coords(1)+i*dx
            end do
        else if (igeometry==2) then
            do i=blk_xlb-1,blk_xub
                r=coords(1)+dx*i
                blk%surf1(i,1,1)=4d0*pi*r*r
                blk%x_interface(i)=coords(1)+i*dx
            end do
            do i=blk_xlb,blk_xub
                r1=coords(1)+dx*(i-1)
                r2=r1+dx
                r=coords(1)+dx*(i-half)
                blk%x_center(i)=r+2*r*dx**2/(12*r**2+dx**2)
                blk%vol(i,1,1)=4d0*pi*(r2**3-r1**3)/3d0
            end do
        end if
    else if (nd==2) then
        coords=blk%pos
        dxyz=blk%dxyz
        call allocate_cell_data_block(blk%vol)
        allocate(blk%x_center(blk_xlb:blk_xub),blk%y_center(blk_ylb:blk_yub),blk%x_interface(blk_xlb-1:blk_xub),blk%y_interface(blk_ylb-1:blk_yub))
        allocate(blk%surf1(blk_xlb-1:blk_xub,blk_ylb:blk_yub,1),blk%surf2(blk_xlb:blk_xub,blk_ylb-1:blk_yub,1))
        allocate(blk%dy(blk_xlb:blk_xub))
        if (igeometry==0) then
            blk%vol=dxyz(1)*dxyz(2)
            blk%surf1=dxyz(2)
            blk%surf2=dxyz(1)
            do i=blk_xlb,blk_xub
                blk%x_center(i)=coords(1)+(i-half)*dxyz(1)
            end do
            do j=blk_ylb,blk_yub
                blk%y_center(j)=coords(2)+(j-half)*dxyz(2)
            end do
            do i=blk_xlb-1,blk_xub
                blk%x_interface(i)=coords(1)+i*dxyz(1)
            end do
            do j=blk_ylb-1,blk_yub
                blk%y_interface(j)=coords(2)+j*dxyz(2)
            end do
        else if (igeometry==1) then
        else if (igeometry==2) then
        end if
    end if
end subroutine initialize_geometry_block

subroutine theta_to_theta(theta,dxyz,theta1,theta2)
    real(8) :: dxyz(3),theta,theta1,theta2
    integer :: i
    i=nint((theta-n_domain(3))/dxyz(2)+half)
    if (i>ny/2) then
        theta=n_domain(3)+dxyz(2)*(ny-i+half)
        theta1=n_domain(3)+dxyz(2)*(ny-i)
        theta2=n_domain(3)+dxyz(2)*(ny-i+1)
    else
        theta=n_domain(3)+dxyz(2)*(i-half)
        theta1=n_domain(3)+dxyz(2)*(i-1)
        theta2=n_domain(3)+dxyz(2)*i
    end if
end subroutine theta_to_theta

subroutine theta_to_theta_interface(theta,dxyz)
    real(8) :: dxyz(3),theta
    integer :: i
    i=nint((theta-n_domain(3))/dxyz(2))
    if (i==ny/2) then
        theta=pi/2d0
    else if (i>ny/2) then
        theta=n_domain(3)+dxyz(2)*(ny-i)
    else
        theta=n_domain(3)+dxyz(2)*i
    end if
end subroutine theta_to_theta_interface

subroutine calculate_sin_cos_theta_interface(theta,sin_val,cos_val,theta_new)
    real(8) :: theta,sin_val,cos_val,dtheta,theta_new
    integer :: n
    dtheta=dxyz(2)
    n=nint((theta-n_domain(3))/dtheta)
    if (n==ny/2) then
        theta_new=pi/2d0
        sin_val=1d0
        cos_val=0d0
    else if (n>ny/2) then
        theta_new=n_domain(3)+dtheta*(ny-n)
        sin_val=sin(theta_new)
        cos_val=-cos(theta_new)
    else
        theta_new=theta
        sin_val=sin(theta)
        cos_val=cos(theta)
    end if
end subroutine calculate_sin_cos_theta_interface

subroutine calculate_sin_cos_theta_center(theta1,theta2,dxyz,sin_val,cos_val,theta)
    !given the bounding theta1 and theta2, 0<theta1<theta2<pi
    !calculate the volume centric sin, cos, store in sin_val and cos_val
    real(8) :: theta1,theta2,dxyz(3),sin_val,cos_val,theta1_temp,theta2_temp,theta_min,theta_max,theta
    logical :: reflect
    if ((theta2-n_domain(3))/dxyz(2)>ny/2) then
        reflect=.true.
    else
        reflect=.false.
    end if
    theta1_temp=theta1
    theta2_temp=theta2
    call theta_to_theta_interface(theta1_temp,dxyz)
    call theta_to_theta_interface(theta2_temp,dxyz)
    theta_min=min(theta1_temp,theta2_temp)
    theta_max=max(theta1_temp,theta2_temp)
    theta=(theta_min*cos(theta_min)-sin(theta_min)-theta_max*cos(theta_max)+sin(theta_max))/(cos(theta_min)-cos(theta_max))
    sin_val=sin(theta)
    if (reflect) then
        cos_val=-cos(theta)
    else
        cos_val=cos(theta)
    end if
end subroutine calculate_sin_cos_theta_center

subroutine calculate_r_center(r,dxyz,r_val)
    !r is the non-weighted center, r_val is the weighted center
    real(8) :: r,dxyz(3),r_val
    r_val=r+(2d0*r*dxyz(1)**2)/(12d0*r**2+dxyz(1)**2)
end subroutine calculate_r_center

function r_to_dr3dr2_ratio(r)
    real(8) :: r,r_to_dr3dr2_ratio
    integer :: i
    i=nint((r-n_domain(1))/dxyz(1)+half)
    r_to_dr3dr2_ratio=dr2dr3_ratio(i)
end function r_to_dr3dr2_ratio

subroutine initialize_hydro_block(blk)
    type(blockdef), pointer :: blk
#if     ischeme==0
    call godunov_initialize_hydro_block(blk)
#elif   ischeme==1
    call godunov_initialize_hydro_block(blk)
#elif   ischeme==2
    call muscl_initialize_hydro_block(blk)
#endif
end subroutine initialize_hydro_block

subroutine godunov_initialize_hydro_block(blk)
    type(blockdef), pointer :: blk
    if (nd==1) then
        call allocate_cell_data_block(blk%w,5)
        call allocate_cell_data_block(blk%u,5)
        call allocate_cell_data_block(blk%dudt_x,5)
        call allocate_cell_data_block(blk%temp)
        call allocate_cell_data_block(blk%egv)
        call allocate_xsurface_data_block(blk%xflux,5)
    else if (nd==2) then
        call allocate_cell_data_block(blk%w,5)
        call allocate_cell_data_block(blk%w1,5)
        call allocate_cell_data_block(blk%w2,5)
        call allocate_cell_data_block(blk%u,5)
        call allocate_cell_data_block(blk%u1,5)
        call allocate_cell_data_block(blk%u2,5)
        call allocate_cell_data_block(blk%dudt_x,5)
        call allocate_cell_data_block(blk%dudt_y,5)
        call allocate_cell_data_block(blk%temp)
        call allocate_cell_data_block(blk%temp1)
        call allocate_cell_data_block(blk%temp2)
        call allocate_cell_data_block(blk%egv)
        call allocate_cell_data_block(blk%egv1)
        call allocate_cell_data_block(blk%egv2)
    end if
end subroutine godunov_initialize_hydro_block

subroutine muscl_initialize_hydro_block(blk)
    type(blockdef), pointer :: blk
    if (nd==1) then
        call allocate_cell_data_block(blk%w,5)
        call allocate_cell_data_block(blk%u,5)
        call allocate_cell_data_block(blk%dudt_x,5)
        call allocate_cell_data_block(blk%temp)
        call allocate_cell_data_block(blk%egv)
        call allocate_cell_data_block(blk%w_xl,5)
        call allocate_cell_data_block(blk%w_xr,5)
        call allocate_cell_data_block(blk%predict_w,5)
        call allocate_cell_data_block(blk%predict_u,5)
        call allocate_cell_data_block(blk%predict_temp)
        call allocate_cell_data_block(blk%predict_egv)
        call allocate_cell_data_block(blk%xslp,5)
        call allocate_xsurface_data_block(blk%xflux,5)
        call allocate_xsurface_data_block(blk%xflux_predict,5)
#if     ieos==2
        call allocate_cell_data_block(blk%cs)
        call allocate_cell_data_block(blk%temp_xl)
        call allocate_cell_data_block(blk%temp_xr)
        call allocate_cell_data_block(blk%egv_xl)
        call allocate_cell_data_block(blk%egv_xr)
#endif
    else if (nd==2) then
        call allocate_cell_data_block(blk%w,5)
        call allocate_cell_data_block(blk%u,5)
        call allocate_cell_data_block(blk%dudt_x,5)
        call allocate_cell_data_block(blk%dudt_y,5)
        call allocate_cell_data_block(blk%temp)
        call allocate_cell_data_block(blk%egv)
        call allocate_cell_data_block(blk%predict_w,5)
        call allocate_cell_data_block(blk%predict_u,5)
        call allocate_cell_data_block(blk%predict_temp)
        call allocate_cell_data_block(blk%predict_egv)
        call allocate_cell_data_block(blk%w0,5)
        call allocate_cell_data_block(blk%u0,5)
        call allocate_cell_data_block(blk%temp0)
        call allocate_cell_data_block(blk%egv0)
        call allocate_cell_data_block(blk%w_xl,5)
        call allocate_cell_data_block(blk%w_xr,5)
        call allocate_cell_data_block(blk%w_yl,5)
        call allocate_cell_data_block(blk%w_yr,5)
        call allocate_cell_data_block(blk%xslp,5)
        call allocate_cell_data_block(blk%yslp,5)
        call allocate_xsurface_data_block(blk%xflux,5)
        call allocate_ysurface_data_block(blk%yflux,5)
#if     ieos==2
        call allocate_cell_data_block(blk%temp_xl)
        call allocate_cell_data_block(blk%temp_xr)
        call allocate_cell_data_block(blk%temp_yl)
        call allocate_cell_data_block(blk%temp_yr)
        call allocate_cell_data_block(blk%egv_xl)
        call allocate_cell_data_block(blk%egv_xr)
        call allocate_cell_data_block(blk%egv_yl)
        call allocate_cell_data_block(blk%egv_yr)
#endif
        if (igeometry==1) then
            call allocate_cell_data_block(blk%am_phi)
            call allocate_cell_data_block(blk%am_phi2)
            call allocate_xsurface_data_block(blk%flux_r,5)
            call allocate_xsurface_data_block(blk%flux_r_amphi)
            call allocate_ysurface_data_block(blk%flux_phi,5)
            call allocate_ysurface_data_block(blk%flux_phi_amphi)
        else if (igeometry==2) then
            call allocate_cell_data_block(blk%vphi)
            call allocate_cell_data_block(blk%vphi2)
            call allocate_cell_data_block(blk%am_theta)
            call allocate_cell_data_block(blk%am_theta2)
            call allocate_cell_data_block(blk%am_phi)
            call allocate_cell_data_block(blk%am_phi2)
            call allocate_cell_data_block(blk%ekphi)
            call allocate_xsurface_data_block(blk%flux_r,5)
            call allocate_xsurface_data_block(blk%flux_r_amphi)
            call allocate_ysurface_data_block(blk%flux_theta,5)
            call allocate_ysurface_data_block(blk%flux_theta_amphi)
        end if
    end if
end subroutine muscl_initialize_hydro_block

subroutine initialize_eos_block(blk)
    type(blockdef), pointer :: blk
#if     ieos==2
    !eos_hllc_analytic only
#if     ieosmodule==1
    call allocate_cell_data_block(blk%H2)
    call allocate_cell_data_block(blk%HI)
    call allocate_cell_data_block(blk%HII)
    call allocate_cell_data_block(blk%electron)
#elif   ieosmodule==2
    call allocate_cell_data_block(blk%HI)
    call allocate_cell_data_block(blk%HII)
#endif
#elif   ieos==3
    !call allocate_cell_data_block(blk%mu)
    !call allocate_cell_data_block(blk%species,7)
#endif
end subroutine initialize_eos_block

subroutine initialize_cooling_block(blk)
    type(blockdef), pointer :: blk
end subroutine initialize_cooling_block

subroutine initialize_radiation_block(blk)
    type(blockdef), pointer :: blk
    if (iradiation==1) then
    else if (iradiation==2) then
        !2d static short characteristics
        call allocate_rays_block(blk%rays,nray_static)
        call allocate_cell_data_block(blk%Erad)
        call allocate_xsurface_data_block(blk%Fradx)
        call allocate_xsurface_data_block(blk%Frady)
        call allocate_cell_data_block(blk%kappa_sc)
        call allocate_cell_data_block(blk%chi_sc)
        call allocate_cell_data_block(blk%rad_source)
    else if (iradiation==3) then
    else if (iradiation==4) then
        if (nd==1) then
            call allocate_xsurface_data_block(blk%kx)
            call allocate_xsurface_data_block(blk%Fradx)
            call allocate_xsurface_data_block(blk%tau)
            call allocate_cell_data_block(blk%cv)
            call allocate_cell_data_block(blk%Erad)
            call allocate_cell_data_block(blk%Erad_int)
            call allocate_cell_data_block(blk%entropy)
            call allocate_cell_data_block(blk%kappa_planck)
            call allocate_cell_data_block(blk%kappa_rosseland)
            call allocate_cell_data_block(blk%sigma_planck)
            call allocate_cell_data_block(blk%sigma_rosseland)
            call allocate_cell_data_block(blk%petsc_egv)
            allocate(blk%nb_coor%xl(1),blk%nb_coor%xu(1))
        else if (nd==2) then
            call allocate_xsurface_data_block(blk%kx)
            call allocate_xsurface_data_block(blk%ky)
            call allocate_xsurface_data_block(blk%Fradx)
            call allocate_xsurface_data_block(blk%Frady)
            call allocate_cell_data_block(blk%cv)
            call allocate_cell_data_block(blk%Erad)
            call allocate_cell_data_block(blk%entropy)
            call allocate_cell_data_block(blk%kappa_planck)
            call allocate_cell_data_block(blk%kappa_rosseland)
            call allocate_cell_data_block(blk%sigma_planck)
            call allocate_cell_data_block(blk%sigma_rosseland)
            allocate(blk%Erad_xl(blk_size_ny),blk%Erad_xu(blk_size_ny),blk%Erad_yl(blk_size_nx),blk%Erad_yu(blk_size_nx),   &
                blk%sigma_rosseland_xl(blk_size_ny),blk%sigma_rosseland_xu(blk_size_ny),    &
                blk%sigma_rosseland_yl(blk_size_nx),blk%sigma_rosseland_yu(blk_size_nx))
            allocate(blk%nb_coor%xl(blk_size_ny),blk%nb_coor%xu(blk_size_ny),blk%nb_coor%yl(blk_size_nx),blk%nb_coor%yu(blk_size_nx))
        end if
    else if (iradiation==5) then
    end if
end subroutine initialize_radiation_block

subroutine initialize_viscosity_block(blk)
    type(blockdef), pointer :: blk
    if (viscous) then
        if (nd==1) then
        else if (nd==2) then
            call allocate_cell_data_block(blk%divv)
            call allocate_xsurface_data_block(blk%vis_tensor_xx)
            call allocate_xsurface_data_block(blk%vis_tensor_xy)
            call allocate_ysurface_data_block(blk%vis_tensor_yy)
            call allocate_ysurface_data_block(blk%vis_tensor_yx)
            if (igeometry==2) then
                call allocate_cell_data_block(blk%vis_tensor_rphi)
                call allocate_cell_data_block(blk%vis_tensor_thetaphi)
                call allocate_cell_data_block(blk%vis_tensor_phitheta)
                call allocate_cell_data_block(blk%vis_tensor_phiphi)
                call allocate_cell_data_block(blk%ek)
                call allocate_cell_data_block(blk%ek2)
                call allocate_cell_data_block(blk%vis_heat)
            end if
        end if
    end if
end subroutine initialize_viscosity_block

subroutine initialize_source_block(blk)
    type(blockdef), pointer :: blk
    real(8) :: r,pos(3),dxyz(3),m
    integer :: i
    if (igeometry==1.or.igeometry==2) then
        call allocate_cell_data_block(blk%s_geo1,5)
        call allocate_cell_data_block(blk%s_geo2,5)
    end if
    if (igravity/=0) then
        if (nd==1) then
            if (igeometry==2) then
                call allocate_cell_data_block(blk%gpotential)
                call allocate_cell_data_block(blk%E_total)
                pos=blk%pos
                dxyz=blk%dxyz
                m=central_star%core%mass
                do i=blk_xlb,blk_xub
                    r=pos(1)+(i-0.5d0)*dxyz(1)
                    blk%gpotential(i,1,1)=-m*gr/r
                end do
            end if
        end if
        call allocate_cell_data_block(blk%s_grav,5)
    end if
    if (igeometry/=0.or.igravity/=0) then
        call allocate_cell_data_block(blk%source,5)
    end if
end subroutine initialize_source_block

subroutine allocate_cell_data_1d(dataarray)
    real(8), allocatable :: dataarray(:,:,:)
    if (nd==1) then
        allocate(dataarray(xlb:xub,1,1))
    else if (nd==2) then
        allocate(dataarray(xlb:xub,ylb:yub,1))
    else
    end if
    dataarray=0d0
end subroutine allocate_cell_data_1d

subroutine allocate_cell_data_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    if (nd==1) then
        allocate(dataarray(xlb:xub,1,1,arraydim))
    else if (nd==2) then
        allocate(dataarray(xlb:xub,ylb:yub,1,arraydim))
    else
    end if
    dataarray=0d0
end subroutine allocate_cell_data_nd

subroutine allocate_cell_data_block_1d(dataarray)
    real(8), allocatable :: dataarray(:,:,:)
    if (nd==1) then
        allocate(dataarray(blk_xlb:blk_xub,1,1))
    else if (nd==2) then
        allocate(dataarray(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    else
    end if
    dataarray=0d0
end subroutine allocate_cell_data_block_1d

subroutine allocate_cell_data_block_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    if (nd==1) then
        allocate(dataarray(blk_xlb:blk_xub,1,1,arraydim))
    else if (nd==2) then
        allocate(dataarray(blk_xlb:blk_xub,blk_ylb:blk_yub,1,arraydim))
    else
    end if
    dataarray=0d0
end subroutine allocate_cell_data_block_nd

subroutine allocate_cell_data_block_front_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    if (nd==1) then
        allocate(dataarray(arraydim,blk_xlb:blk_xub,1,1))
    else if (nd==2) then
        allocate(dataarray(arraydim,blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    else
    end if
    dataarray=0d0
end subroutine allocate_cell_data_block_front_nd

subroutine allocate_xsurface_data_1d(dataarray)
    real(8), allocatable :: dataarray(:,:,:)
    if (nd==1) then
        allocate(dataarray(xlb:xub-1,1,1))
    else if (nd==2) then
        allocate(dataarray(xlb:xub-1,ylb:yub,1))
    else
    end if
    dataarray=0d0
end subroutine allocate_xsurface_data_1d

subroutine allocate_xsurface_data_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    if (nd==1) then
        allocate(dataarray(xlb:xub-1,1,1,arraydim))
    else if (nd==2) then
        allocate(dataarray(xlb:xub-1,ylb:yub,1,arraydim))
    end if
    dataarray=0d0
end subroutine allocate_xsurface_data_nd

subroutine allocate_xsurface_data_block_1d(dataarray)
    real(8), allocatable :: dataarray(:,:,:)
    if (nd==1) then
        allocate(dataarray(blk_xlb:blk_xub-1,1,1))
    else if (nd==2) then
        allocate(dataarray(blk_xlb:blk_xub-1,blk_ylb:blk_yub,1))
    else
    end if
    dataarray=0d0
end subroutine allocate_xsurface_data_block_1d

subroutine allocate_xsurface_data_block_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    if (nd==1) then
        allocate(dataarray(blk_xlb:blk_xub-1,1,1,arraydim))
    else
        allocate(dataarray(blk_xlb:blk_xub-1,blk_ylb:blk_yub,1,arraydim))
    end if
    dataarray=0d0
end subroutine allocate_xsurface_data_block_nd

subroutine allocate_xsurface_data_block_front_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    if (nd==1) then
        allocate(dataarray(arraydim,blk_xlb:blk_xub-1,1,1))
    else
        allocate(dataarray(arraydim,blk_xlb:blk_xub-1,blk_ylb:blk_yub,1))
    end if
    dataarray=0d0
end subroutine allocate_xsurface_data_block_front_nd

subroutine allocate_ysurface_data_1d(dataarray)
    real(8), allocatable :: dataarray(:,:,:)
    if (nd==1) then
        !1d is x direction only
    else if (nd==2) then
        allocate(dataarray(xlb:xub,ylb:yub-1,1))
    else
    end if
    dataarray=0d0
end subroutine allocate_ysurface_data_1d

subroutine allocate_ysurface_data_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    if (nd==1) then
        !1d is x direction only
    else if (nd==2) then
        allocate(dataarray(xlb:xub,ylb:yub-1,1,arraydim))
    else
    end if
    dataarray=0d0
end subroutine allocate_ysurface_data_nd

subroutine allocate_ysurface_data_block_1d(dataarray)
    real(8), allocatable :: dataarray(:,:,:)
    if (nd==1) then
        !1d is x direction only
    else if (nd==2) then
        allocate(dataarray(blk_xlb:blk_xub,blk_ylb:blk_yub-1,1))
    else
    end if
    dataarray=0d0
end subroutine allocate_ysurface_data_block_1d

subroutine allocate_ysurface_data_block_nd(dataarray,arraydim)
    real(8), allocatable :: dataarray(:,:,:,:)
    integer :: arraydim
    if (nd==1) then
        !1d is x direction only
    else if (nd==2) then
        allocate(dataarray(blk_xlb:blk_xub,blk_ylb:blk_yub-1,1,arraydim))
    else
    end if
    dataarray=0d0
end subroutine allocate_ysurface_data_block_nd

subroutine allocate_rays_block(rays,nray)
    !each cell center has nray
    type(ray2ddef), dimension(:,:,:), allocatable :: rays
    integer :: nray
    allocate(rays(nray,blk_size_nx,blk_size_ny))
end subroutine allocate_rays_block

subroutine initialize_global_parameters_int(pvalue,p,printout)
    !give pvalue to p, p is protected
    !only this subroutine can pass values to global parameters
    integer :: pvalue,p
    character(len=*) :: printout
    p=pvalue
    if (rank==0) write(*,'(A32,I18)') printout,p
end subroutine initialize_global_parameters_int

subroutine initialize_global_parameters_real(pvalue,p,printout)
    !give pvalue to p, p is protected
    !only this subroutine can pass values to global parameters
    real(8) :: pvalue,p
    character(len=*) :: printout
    p=pvalue
    if (rank==0) write(*,'(A32,ES18.8E3)') printout,p
end subroutine initialize_global_parameters_real

subroutine initialize_global_parameters_logical(pvalue,p,printout)
    !give pvalue to p, p is protected
    !only this subroutine can pass values to global parameters
    logical :: pvalue,p
    character(len=*) :: printout
    p=pvalue
    if (rank==0) write(*,'(A32,L18)') printout,p
end subroutine initialize_global_parameters_logical

subroutine initialize_mesh()
#if     ischeme==0
    n_hydro_guard=2
#elif   ischeme==1
    n_hydro_guard=1
#elif   ischeme==2
    n_hydro_guard=2
#elif   ischeme==3
    n_hydro_guard=2
#endif
    if (nd==1) then
        xlb=1-n_hydro_guard
        xub=nx+n_hydro_guard
        ylb=1
        yub=1
    else if (nd==2) then
        xlb=1-n_hydro_guard
        xub=nx+n_hydro_guard
        ylb=1-n_hydro_guard
        yub=ny+n_hydro_guard
    end if
    if (rank==0) then
        write(*,'(A20,I15)') 'n_hydro_guard=',n_hydro_guard
    end if
end subroutine initialize_mesh

subroutine blk_traversal(sub)
    type(blockdef), pointer :: blk
    procedure(blk_operator) :: sub
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        call sub(blk)
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end subroutine blk_traversal

subroutine ijk_to_coords(ijk,blk_dxyz,blk_coords,coords)
    !convert ijk to coordinates
    integer :: ijk(3)
    real(8), dimension(3) :: coords,blk_dxyz,blk_coords
    real(8) :: theta1,theta2,r_val,sin_val,cos_val
    if (nd==1) then
        coords(1)=blk_coords(1)+(ijk(1)-half)*blk_dxyz(1)
        coords(2)=0d0
        coords(3)=0d0
    else if (nd==2) then
        coords(1)=blk_coords(1)+(ijk(1)-half)*blk_dxyz(1)
        coords(2)=blk_coords(2)+(ijk(2)-half)*blk_dxyz(2)
        coords(3)=0d0
        if (igeometry==2) then
            call theta_to_theta(coords(2),blk_dxyz,theta1,theta2)
            call calculate_sin_cos_theta_center(theta1,theta2,dxyz,sin_val,cos_val,coords(2))
            call calculate_r_center(coords(1),dxyz,r_val)
            coords(1)=r_val
        end if
    else
    end if
end subroutine ijk_to_coords

subroutine coords_to_ijk(coords,blk_dxyz,blk_coords,ijk)
    !converts coordnates to ijk
    integer :: ijk(3)
    real(8), dimension(3) :: coords,blk_dxyz,blk_coords
    if (nd==1) then
        ijk(1)=ceiling((coords(1)-blk_coords(1))/blk_dxyz(1))
        ijk(2)=1
        ijk(3)=1
    else if (nd==2) then
        ijk(1)=ceiling((coords(1)-blk_coords(1))/blk_dxyz(1))
        ijk(2)=ceiling((coords(2)-blk_coords(2))/blk_dxyz(2))
        ijk(3)=1
    end if
end subroutine coords_to_ijk

subroutine collective_sub(sub1,sub2,sub3,outcome)
    !the subroutine that carries out the collective operation
    !sub1 operates on each block of every processor
    !sub2 operates on each processor, save the result to processor%collective_processor_result
    !sub3 operates on rank0
    type(blockdef), pointer :: blk
    procedure(collective_oper_sub1) :: sub1
    procedure(collective_oper_sub2) :: sub2
    procedure(collective_oper_sub3) :: sub3
    real(8) :: blk_result,outcome
    real(8), dimension(:), allocatable :: blk_result_array
    integer :: i
    processor%collective_processor_result=0d0
    processor%collective_array=0d0
    processor%collective_result=0d0
    allocate(blk_result_array(np_nblk(rank+1)))
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        call sub1(blk,blk_result)
        blk_result_array(i)=blk_result
        blk=>blk%blk_next
    end do
    call sub2(blk_result_array)
    call sub3()
    outcome=processor%collective_result
    deallocate(blk_result_array)
end subroutine collective_sub

subroutine collective_processor_min(blk_result_array)
    real(8), dimension(:), allocatable :: blk_result_array
    processor%collective_processor_result=minval(blk_result_array)
end subroutine collective_processor_min

subroutine collective_processor_sum(blk_result_array)
    real(8), dimension(:), allocatable :: blk_result_array
    processor%collective_processor_result=sum(blk_result_array)
end subroutine collective_processor_sum

subroutine collective_rank0_min()
    integer :: ierr
    call mpi_gather(processor%collective_processor_result,1,MPI_REAL8,processor%collective_array,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    processor%collective_result=minval(processor%collective_array)
    call mpi_bcast(processor%collective_result,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine collective_rank0_min

subroutine collective_rank0_sum()
    integer :: ierr
    call mpi_gather(processor%collective_processor_result,1,MPI_REAL8,processor%collective_array,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    processor%collective_result=sum(processor%collective_array)
    call mpi_bcast(processor%collective_result,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine collective_rank0_sum

subroutine check_global_parameter_consistency()
    integer :: max_static_level
    character(len=256) :: alert
#if     ieos==2
    if (opacity_gas_rho_min<=environment%rhomin) then
        alert='ieos=2 and opacity_gas_rho_min<=environment%rhomin'
        call abort_achilles(alert)
    end if
    if (opacity_gas_rho_max>=environment%rhomax) then
        alert='ieos=2 and opacity_gas_rho_max>=environment%rhomax'
        call abort_achilles(alert)
    end if
#endif
    if (refine_type=='static'.or.refine_type=='mixed') then
        max_static_level=maxval(refine_region(:,5))
        if (max_refine_level<max_static_level) then
            alert='max_refine_level<max_static_level'
            call abort_achilles(alert)
        end if
    end if
    if (nd/=1) then
        if (blk_size_nx/=blk_size_ny) then
            alert='block must be square'
            call abort_achilles(alert)
        end if
        if (hydro_bound_type(1)==3.or.hydro_bound_type(2)==3) then
            if (hydro_bound_type(1)/=3.or.hydro_bound_type(2)/=3) then
                alert='xl and xu must be periodic at the same time'
                call abort_achilles(alert)
            end if
        end if
        if (hydro_bound_type(3)==3.or.hydro_bound_type(4)==3) then
            if (hydro_bound_type(3)/=3.or.hydro_bound_type(4)/=3) then
                alert='yl and yu must be periodic at the same time'
                call abort_achilles(alert)
            end if
        end if
        if (igeometry==1) then
            alert='using cylindrical coordinate'
            call print_achilles(alert)
            if (n_domain(1)<=0d0) then
                alert='r inner domain <= 0'
                call abort_achilles(alert)
            end if
            if (n_domain(3)<0d0) then
                alert='phi lower domain < 0.'
                call abort_achilles(alert)
            end if
            if (n_domain(4)>2d0*pi) then
                alert='phi upper domain > 2pi.'
                call abort_achilles(alert)
            end if
        else if (igeometry==2) then
            alert='using polar coordinate'
            call print_achilles(alert)
            if (n_domain(1)<=0d0) then
                alert='r inner domain <= 0'
                call abort_achilles(alert)
            end if
            if (n_domain(3)<=0d0) then
                alert='theta lower domain <= 0.'
                call abort_achilles(alert)
            end if
            if (n_domain(4)>=pi) then
                alert='theta upper domain >= pi.'
                call abort_achilles(alert)
            end if
        end if
    end if
end subroutine check_global_parameter_consistency

subroutine diagnostics()
end subroutine diagnostics

subroutine check_nan()
    type(blockdef), pointer :: blk
    integer :: i,j
    character(len=128) :: str
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        do j=blk_xlb,blk_xub
            if (isnan(blk%w(j,1,1,1)).or.blk%w(j,1,1,1)<=0) then
                print *,'density',time_sys%ntimestep,rank,i,j,blk%key,blk%level
                str='nan or non-positive'
                call abort_achilles(str)
            end if
            if (isnan(blk%temp(j,1,1))) then
                print *,'temp',time_sys%ntimestep,rank,i,j,blk%key,blk%level
                str='nan'
                call abort_achilles(str)
            end if
            if (isnan(blk%egv(j,1,1))) then
                print *,'egv',time_sys%ntimestep,rank,i,j,blk%key,blk%level
                str='nan'
                call abort_achilles(str)
            end if
            if (blk%w(j,1,1,5)<0) then
                print *,'pres',time_sys%ntimestep,rank,i,j,blk%blk_id,blk%blk_next%blk_id
                print *,blk%w(j,1,1,5)
                print *,blk%blk_next%w(:,1,1,5)
                print *,blk%level,blk%blk_next%level
                str='negative'
                call abort_achilles(str)
            end if
        end do
        blk=>blk%blk_next
    end do
end subroutine check_nan

subroutine abort_achilles(str)
    integer :: ierr
    character(len=128), intent(in) :: str
    if (rank==0) then
        print *,trim(str)
    end if
    call mpi_finalize(ierr)
    stop
end subroutine abort_achilles

subroutine print_achilles(str)
    integer :: ierr
    character(len=128), intent(in) :: str
    if (rank==0) then
        print *,trim(str)
    end if
end subroutine print_achilles

end module datastructure
