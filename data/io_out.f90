module io_out
use datastructure
use communication
use mathlib
use phylib
use eos
use HDF5
use mpi
implicit none

character(len=32) :: dsetname(100),filename_head
integer :: output_var_number

interface output_simulation_attr
    module procedure output_simulation_attr_double,output_simulation_attr_int
end interface output_simulation_attr

contains

subroutine simulation_config()
    character(len=24) :: fmt1,fmt2,fmt3
    fmt1='(A20,1I15)'
    fmt2='(A20,1ES15.5E2)'
    fmt3='(A20,L15)'
    if (rank==0) then
        call chdir(path_out)
        open(unit=14,file='configuration.dat',status='replace',action='write')
        call print_config(14)
        close(14)
        call chdir(path_root)
    end if
end subroutine simulation_config

subroutine print_config(i)
    integer :: i
    character(len=24) :: fmt1,fmt2,fmt3
    fmt1='(A20,1I15)'
    fmt2='(A20,1ES15.5E2)'
    fmt3='(A20,L15)'
    write(unit=i,fmt=*) '&config_nml'
    write(unit=i,fmt=fmt2) 'CFL_number=', CFL
    write(unit=i,fmt=fmt1) 'nd=', nd
    write(unit=i,fmt=fmt1) 'number_of_frame=', nframe
    write(unit=i,fmt=fmt2) 'coord1min=', n_domain(1)
    write(unit=i,fmt=fmt2) 'coord1max=', n_domain(2)
    write(unit=i,fmt=fmt2) 'coord2min=', n_domain(3)
    write(unit=i,fmt=fmt2) 'coord2max=', n_domain(4)
    write(unit=i,fmt=fmt1) 'ncoord1=', nx
    write(unit=i,fmt=fmt1) 'ncoord2=', ny
    write(unit=i,fmt=fmt1) 'ncoord3=', nz
    write(unit=i,fmt=fmt1) 'igravity=', igravity
    write(unit=i,fmt=fmt1) 'igeometry=', igeometry
    write(unit=i,fmt=fmt1) 'icooling=', icooling
    write(unit=i,fmt=fmt1) 'iradiation=', iradiation
    write(unit=i,fmt=fmt1) 'maximum logical level', llevel_max
    write(unit=i,fmt=fmt2) 'mean atomic weight=', maw
    write(unit=i,fmt=fmt2) 'gamma_gas=', gamma_gas 
    write(unit=i,fmt=*) '/'
end subroutine print_config

subroutine output_blocks()
    !rank0 create the output file and write the meta data to the file
    !each rank output its own blocks
    character(len=32) :: filename_hdf,filename_xml
    character(len=16) :: fmt1,s1
    integer :: i,ierr
    fmt1='(I5.5)'
    write(s1,fmt1) iframe
    filename_hdf=trim(filename_head)//trim(s1)//'.h5'
    filename_xml=trim(filename_head)//trim(s1)//'.xmf'
    call calculate_output_data()
    if (rank==0) call output_meta_data(filename_hdf)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    do i=0,np-1
        call hdf_output_block_sharedfile(filename_hdf,i)
        call mpi_barrier(MPI_COMM_WORLD,ierr)
    end do
    if (nd==2) then
        if (rank==0) then
            call create_xdmf_xml(filename_hdf,filename_xml)
        end if
    end if
end subroutine output_blocks

subroutine calculate_output_data()
    !calculate some derivative data here
    type(blockdef), pointer :: blk
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        if (iradiation/=0) then
            call calculate_entropy_block(blk)
        end if
        blk=>blk%blk_next
    end do
end subroutine calculate_output_data

subroutine output_meta_data(filename)
    type(blockdef), pointer :: blk_traversal
    integer :: i,j,ijk(3),error
    character(len=16) :: fmt1,s1
    character(len=32) :: filename
    character(len=128) :: alert,att_str
    integer(HID_T) :: file_id
    call chdir(path_out)
    call h5open_f(error)
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)
    call output_simulation_attr(file_id,'time',time_sys%t)
    call output_simulation_attr(file_id,'nblocks',nblk_total)
    call output_simulation_attr(file_id,'blk_size_nx',blk_size_nx)
    call save_tree(file_id)
    call h5fclose_f(file_id,error)
    call h5close_f(error)
    call chdir(path_root)
end subroutine output_meta_data

subroutine dsetname_filename()
    character(len=32) :: filename2
    integer :: i
    call chdir(path_out)
    open(unit=14,file='../output_var_info.dat',status='old',action='read')
    read(unit=14,fmt=*) output_var_number
    do i=1,output_var_number
        read(unit=14,fmt=*) dsetname(i)
    end do
    if (igeometry==2) then
        if (nd==2) then
            dsetname(output_var_number+1)='vphi'
            output_var_number=output_var_number+1
        end if
    end if
    if (iradiation==1) then
    else if (iradiation==3) then
    else if (iradiation==4) then
        if (nd==1) then
            dsetname(output_var_number+1)='Erad'
            dsetname(output_var_number+2)='Fradx'
            dsetname(output_var_number+3)='Rosseland'
            dsetname(output_var_number+4)='Planck'
            dsetname(output_var_number+5)='entropy'
            output_var_number=output_var_number+5
        else if (nd==2) then
            dsetname(output_var_number+1)='Erad'
            dsetname(output_var_number+2)='Fradx'
            dsetname(output_var_number+3)='sigma_Rosseland'
            dsetname(output_var_number+4)='sigma_Planck'
            output_var_number=output_var_number+4
        end if
    else if (iradiation==5) then
        !if (nd==1) then
        !    dsetname(output_var_number+1)='Erad_mg'
        !    dsetname(output_var_number+2)='Rosseland_mg'
        !    dsetname(output_var_number+3)='Planck_mg'
        !    dsetname(output_var_number+4)='entropy'
        !    output_var_number=output_var_number+4
        !end if
    end if
#if     ieos==2
#if     ieosmodule==1
    dsetname(output_var_number+1)='H2'
    dsetname(output_var_number+2)='HI'
    dsetname(output_var_number+3)='HII'
    output_var_number=output_var_number+3
#elif   ieosmodule==2
    dsetname(output_var_number+1)='HI'
    dsetname(output_var_number+2)='HII'
    output_var_number=output_var_number+2
#endif
#endif
    read(unit=14,fmt=*) filename2
    close(14)
    filename_head=trim(filename2)
    call chdir(path_root)
end subroutine dsetname_filename

subroutine hdf_output_block_sharedfile(filename,irank)
    !hdf output supports 1d and 2d
    !each blocks is a group in the hdf file
    type(blockdef), pointer :: blk_traversal
    integer :: i,j,k,ijk(3)
    integer :: error,dset_rank,arank,irank
    integer(HSIZE_T), dimension(:), allocatable :: dims,dims_xsurface,dims_ysurface,dims_mg
    integer(HSIZE_T) :: blk_xdim(1),blk_ydim(1),att_dim(1)
    real(8) :: rho,temp,p,p_temp,log10_temp,log10_egv,eg,cs,xyz(3),rthetaphi(3),egv,nhtot,species4(4)
    real(8), allocatable :: species5(:)
    real(8), dimension(:,:), allocatable :: cell_data,x_surfdata,y_surfdata
    real(8), dimension(:,:,:), allocatable :: cell_data_mg,x_surfdata_mg
    character(len=16) :: fmt1,s1
    character(len=32) :: filename,groupname,blockname,attname
    character(len=128) :: alert,att_str
    integer(HID_T) :: file_id,group_id,dset_id,dspace_id,aspace_id,attr_id
    if (rank==irank) then
        call chdir(path_out)
        call h5open_f(error)
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
        dset_rank=2
        allocate(dims(dset_rank),dims_xsurface(dset_rank))
        dims=(/blk_size_nx,blk_size_ny/)
        allocate(cell_data(dims(1),dims(2)))
        if (iradiation==5) then
            !allocate(dims_mg(dset_rank+1))
            !dims_mg=(/nmg,blk_size_nx,blk_size_ny/)
            !allocate(cell_data_mg(dims_mg(1),dims_mg(2),dims_mg(3)))
        end if
        blk_xdim=blk_size_nx+1
        dims_xsurface(1)=blk_size_nx+1
        if (nd==1) then
            blk_ydim=1
            dims_xsurface(2)=1
        else if (nd==2) then
            blk_ydim=blk_size_ny+1
            dims_xsurface(2)=blk_size_ny
        end if
        allocate(x_surfdata(dims_xsurface(1),dims_xsurface(2)))
        blk_traversal=>blk_processor_head
        do k=1,np_nblk(rank+1)
#if     ieos==2
            !analytic realistic gas EOS
            if (nd==1) then
#if     ieosmodule==1
                allocate(species5(5))
                do i=1,blk_size_nx
                    rho=blk_traversal%w(i,1,1,1)
                    temp=blk_traversal%temp(i,1,1)
                    call calspecies(rho,temp,species5)
                    nhtot=2*species5(1)+species5(2)+species5(3)
                    blk_traversal%H2(i,1,1)=2*species5(1)/nhtot
                    blk_traversal%HI(i,1,1)=species5(2)/nhtot
                    blk_traversal%HII(i,1,1)=species5(3)/nhtot
                    blk_traversal%cs(i,1,1)=calculate_adiabatic_cs(rho,temp)
                end do
                deallocate(species5)
#elif   ieosmodule==2
                do i=1,blk_size_nx
                    rho=blk_traversal%w(i,1,1,1)
                    temp=blk_traversal%temp(i,1,1)
                    call calspecies(rho,temp,species4)
                    nhtot=species4(1)+species4(2)
                    blk_traversal%HI(i,1,1)=species4(1)/nhtot
                    blk_traversal%HII(i,1,1)=species4(2)/nhtot
                    blk_traversal%cs(i,1,1)=calculate_adiabatic_cs(rho,temp)
                end do
#endif
            else if (nd==2) then
            end if
#endif
            fmt1='(I5.5)'
            write(s1,fmt1) blk_traversal%blk_id
            groupname='blk'//trim(s1)
            call h5gcreate_f(file_id,groupname,group_id,error)
            !create the datasets
            if (nd==1) then
                call h5screate_simple_f(1,blk_xdim,dspace_id,error)
                call h5dcreate_f(group_id,'mesh_x',h5t_native_double,dspace_id,dset_id,error)
                arank=1
                att_dim(1)=1
                call h5screate_simple_f(arank,att_dim,aspace_id,error)
                call h5acreate_f(dset_id,'level',h5t_native_integer,aspace_id,attr_id,error)
                call h5awrite_f(attr_id,h5t_native_integer,blk_traversal%level,att_dim,error)
                call h5aclose_f(attr_id,error)
                call h5sclose_f(aspace_id, error)
                call h5dwrite_f(dset_id,h5t_native_double,blk_traversal%mesh_x,blk_xdim,error)
                call h5dclose_f(dset_id,error)
            else if (nd==2) then
                call h5screate_simple_f(1,blk_xdim,dspace_id,error)
                call h5dcreate_f(group_id,'mesh_x',h5t_native_double,dspace_id,dset_id,error)
                call h5dwrite_f(dset_id,h5t_native_double,blk_traversal%mesh_x,blk_xdim,error)
                call h5dclose_f(dset_id,error)
                call h5screate_simple_f(1,blk_ydim,dspace_id,error)
                call h5dcreate_f(group_id,'mesh_y',h5t_native_double,dspace_id,dset_id,error)
                call h5dwrite_f(dset_id,h5t_native_double,blk_traversal%mesh_y,blk_ydim,error)
                call h5dclose_f(dset_id,error)
            end if
            do i=1,output_var_number
                select case(dsetname(i))
                !some variables are cell centered quantities, some are surface quantities or a data array for each cell center.
                case('rho')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%w(1:blk_size_nx,1:blk_size_ny,1,irho)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('vx')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%w(1:blk_size_nx,1:blk_size_ny,1,ivx)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('vy')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%w(1:blk_size_nx,1:blk_size_ny,1,ivy)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('vz')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%w(1:blk_size_nx,1:blk_size_ny,1,ivz)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('pres')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%w(1:blk_size_nx,1:blk_size_ny,1,ipres)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('cs')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%cs(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('temp')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%temp(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('egv')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%egv(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('entropy')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%entropy(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('vphi')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%vphi(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('H2')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%H2(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('HI')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%HI(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('HII')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%HII(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('Rosseland')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%kappa_rosseland(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('sigma_Rosseland')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%sigma_rosseland(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('Rosseland_mg')
                    call h5screate_simple_f(dset_rank+1,dims_mg,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data_mg=blk_traversal%kappa_rosseland_mg(1:nmg,:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data_mg,dims_mg,error)
                case('Planck')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%kappa_planck(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('sigma_Planck')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data=blk_traversal%sigma_planck(1:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('Planck_mg')
                    call h5screate_simple_f(dset_rank+1,dims_mg,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data_mg=blk_traversal%kappa_planck_mg(1:nmg,:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data_mg,dims_mg,error)
                case('Erad')
                    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    if (resolve_rad_dt) then
                        cell_data=blk_traversal%Erad_int(1:blk_size_nx,1:blk_size_ny,1)
                    else
                        cell_data=blk_traversal%Erad(1:blk_size_nx,1:blk_size_ny,1)
                    end if
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data,dims,error)
                case('Erad_mg')
                    call h5screate_simple_f(dset_rank+1,dims_mg,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    cell_data_mg=blk_traversal%Erad_mg(1:nmg,:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,cell_data_mg,dims_mg,error)
                case('Fradx')
                    call h5screate_simple_f(dset_rank,dims_xsurface,dspace_id,error)
                    call h5dcreate_f(group_id,trim(dsetname(i)),h5t_native_double,dspace_id,dset_id,error)
                    x_surfdata=blk_traversal%Fradx(0:blk_size_nx,1:blk_size_ny,1)
                    call h5dwrite_f(dset_id,h5t_native_double,x_surfdata,dims_xsurface,error)
                case default
                    alert='no such dsetname:'
                    print *,trim(groupname)//trim(dsetname(i))
                    call abort_achilles(alert)
                end select
                call h5dclose_f(dset_id,error)
                call h5sclose_f(dspace_id,error)
            end do
            CALL h5gclose_f(group_id,error)
            blk_traversal=>blk_traversal%blk_next
        end do
        deallocate(dims,dims_xsurface)
        deallocate(cell_data,x_surfdata)
        if (iradiation==5) then
            !deallocate(dims_mg,cell_data_mg)
        end if
        call h5fclose_f(file_id,error)
        call h5close_f(error)
        call chdir(path_root)
    end if
end subroutine hdf_output_block_sharedfile

subroutine output_simulation_attr_double(file_id,str,val)
    integer(HID_T) :: file_id,aspace_id,attr_id
    integer(HSIZE_T) :: att_dim(1)
    integer :: arank,error
    real(8) :: val
    character(len=*) :: str
    arank=1
    att_dim(1)=1
    call h5screate_simple_f(arank,att_dim,aspace_id,error)
    call h5acreate_f(file_id,trim(str),h5t_native_double,aspace_id,attr_id,error)
    call h5awrite_f(attr_id,h5t_native_double,val,att_dim,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id, error)
end subroutine output_simulation_attr_double

subroutine output_simulation_attr_int(file_id,str,val)
    integer(HID_T) :: file_id,aspace_id,attr_id
    integer(HSIZE_T) :: att_dim(1)
    integer :: arank,error,val
    character(len=*) :: str
    arank=1
    att_dim(1)=1
    call h5screate_simple_f(arank,att_dim,aspace_id,error)
    call h5acreate_f(file_id,trim(str),h5t_native_integer,aspace_id,attr_id,error)
    call h5awrite_f(attr_id,h5t_native_integer,val,att_dim,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id, error)
end subroutine output_simulation_attr_int

subroutine save_tree(file_id)
    integer(HID_T) :: file_id,dspace_id,dset_id
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer, dimension(:,:), allocatable :: keys
    integer :: dset_rank,error,i
    dset_rank=2
    call count_nblocks_tree()
    call output_simulation_attr(file_id,'nblocks_tree',nblocks_tree)
    allocate(dims(2),keys(nblocks_tree,6))
    keys=0
    dims=(/nblocks_tree,6/)
    call save_the_keys(keys)
    call h5screate_simple_f(dset_rank,dims,dspace_id,error)
    call h5dcreate_f(file_id,'block_tree',h5t_native_integer,dspace_id,dset_id,error)
    call h5dwrite_f(dset_id,h5t_native_integer,keys,dims,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)
    deallocate(dims,keys)
end subroutine save_tree

subroutine count_nblocks_tree()
    integer :: n
    n=1
    call count_child(blk_root,n)
    nblocks_tree=n
end subroutine count_nblocks_tree

subroutine count_child(blk,n)
    !preorder traversal
    type(blockdef), pointer :: blk
    integer :: n
    if (nd==1) then
        if (associated(blk%blk_xl)) then
            n=n+1
            call count_child(blk%blk_xl,n)
        end if
        if (associated(blk%blk_xu)) then
            n=n+1
            call count_child(blk%blk_xu,n)
        end if
    else if (nd==2) then
    end if
end subroutine count_child

subroutine save_the_keys(keys)
    integer, dimension(:,:), allocatable :: keys
    integer :: i
    i=1
    call save_key(blk_root,keys,i)
end subroutine save_the_keys

subroutine save_key(blk,keys,i)
    !preorder traversal
    type(blockdef), pointer :: blk
    integer, dimension(:,:), allocatable :: keys
    integer :: i,key_p(3),key(3)
    if (associated(blk%blk_p)) then
        key_p=blk%blk_p%key
    else
        key_p=blk%key
    end if
    key=blk%key
    keys(i,1:6)=(/key_p(1:3),key(1:3)/)
    if (nd==1) then
        if (associated(blk%blk_xl)) then
            i=i+1
            call save_key(blk%blk_xl,keys,i)
        end if
        if (associated(blk%blk_xu)) then
            i=i+1
            call save_key(blk%blk_xu,keys,i)
        end if
    else if (nd==2) then
    end if
end subroutine save_key

subroutine create_xdmf_xml(filename_hdf,filename_xml)
    character(len=32) :: fmt1,s1,filename_hdf,filename_xml,blockname
    type(blockdef), pointer :: blk_traversal
    integer :: i,j,k
    call chdir(path_out)
    fmt1='(I5.5)'
    open(unit=15,file=filename_xml,status='replace',action='write')
    write(15,'(a)')'<Xdmf Version="2.0">'
    write(15,'(a)')'<Domain Name="domain">'
    write(15,'(a)')'<Grid Name="mesh" GridType="Collection">'
    blk_traversal=>blk_head
    do while (associated(blk_traversal))
        write(s1,fmt1) blk_traversal%blk_id
        blockname='blk'//trim(s1)
        write(15,*)'    <Grid Name="',trim(blockname),'" GridType="Uniform">'
        write(15,*)'    <Topology TopologyType="2DRectMesh" NumberOfElements="',blk_size_nx+1,blk_size_ny+1,'"/>'
        write(15,*)'    <Geometry GeometryType="VXVY">'
        write(15,*)'        <DataItem Dimensions="',blk_size_nx+1,'" NumberType="Double" Precision="8" Format="HDF">'
        write(15,*)'            '//trim(filename_hdf)//':/'//trim(blockname)//'/'//'mesh_x'
        write(15,*)'        </DataItem>'
        write(15,*)'        <DataItem Dimensions="',blk_size_ny+1,'" NumberType="Double" Precision="8" Format="HDF">'
        write(15,*)'            '//trim(filename_hdf)//':/'//trim(blockname)//'/'//'mesh_y'
        write(15,*)'        </DataItem>'
        write(15,*)'    </Geometry>'
        do i=1,output_var_number
            write(15,*)'    <Attribute Name="',trim(dsetname(i)),'" AttributeType="Scalar" Center="Cell">'
            write(15,*)'        <DataItem Dimensions="',blk_size_nx,blk_size_ny,'" NumberType="Double" Precision="8" Format="HDF">'
            write(15,*)'        '//trim(filename_hdf)//':/'//trim(blockname)//'/'//trim(dsetname(i))
            write(15,*)'        </DataItem>'
            write(15,*)'    </Attribute>'
        end do
        write(15,*)'    </Grid>'
        blk_traversal=>blk_traversal%blk_next
    end do
    write(15,'(a)')'</Grid>'
    write(15,'(a)')'</Domain>'
    write(15,'(a)')'</Xdmf>'
    close(15)
    call chdir(path_root)
end subroutine create_xdmf_xml

subroutine record(record_array,record_i)
    !append record_array to record_history
    real(8), allocatable :: record_array(:),record_history_swap(:,:)
    integer :: record_i,record_chunk_size
    if (rank==0) then
        record_chunk_size=1000
        record_array_size=size(record_array)
        if (allocated(record_history)) then
            if (mod(record_i,record_chunk_size)==1) then    !overfill, need to expand record_history
                record_chunk_i=record_chunk_i+1
                allocate(record_history_swap(record_chunk_i*record_chunk_size,record_array_size))
                record_history_swap=0d0
                record_history_swap(1:(record_chunk_i-1)*record_chunk_size,:)=record_history
                deallocate(record_history)
                allocate(record_history(record_chunk_i*record_chunk_size,record_array_size))
                record_history=0d0
                record_history=record_history_swap
                record_history(record_i,:)=record_array
                write(*,fmt='(32ES14.3E3)') record_history(record_i,:)
                deallocate(record_history_swap)
            else                                            !old record_history can still uphold new data
                record_history(record_i,:)=record_array
                write(*,fmt='(32ES14.3E3)') record_history(record_i,:)
            end if
        else   !the first record
            record_i=1
            record_chunk_i=1
            record_array_size=size(record_array)
            allocate(record_history(record_chunk_i*record_chunk_size,record_array_size))
            record_history(record_i,:)=record_array
            write(*,fmt='(32ES14.3E3)') record_history(record_i,:)
        end if
    end if
end subroutine record

subroutine save_recorded_data(nrecord)
    integer :: i,nrecord,length
    character(len=32) :: fmt1,s1
    if (rank==0) then
        length=size(record_history,2)
        call chdir(path_out)
        open(unit=16,file='history.data',status='replace',action='write')
        do i=1,nrecord
            write(unit=16,fmt='(32ES15.3E3)') record_history(i,:)
        end do
        close(16)
        call chdir(path_root)
    end if
end subroutine save_recorded_data

end module io_out
