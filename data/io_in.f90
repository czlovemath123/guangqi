module io_in
use datastructure
use communication
use io_out
use mathlib
use phylib
use boundary
use eos
use problem
use HDF5
use mpi
implicit none

contains

subroutine hdf_read()
    character(len=32) :: filename_hdf,blockname,attname
    integer(HID_T) :: file_id,att_id
    integer(HSIZE_T) :: att_dim(1)
    integer :: error,ierr
    character(len=16) :: fmt1,s1
    if (rank==0) then
        fmt1='(I5.5)'
        write(s1,fmt1) restart_iframe
        filename_hdf=trim(filename_head)//trim(s1)//'.h5'
        call chdir(path_out)
        att_dim=1
        call h5open_f(error)
        call h5fopen_f(filename_hdf,H5F_ACC_RDWR_F,file_id,error)
        call h5aopen_f(file_id,'nblocks_tree',att_id,error)
        call h5aread_f(att_id,h5t_native_integer,nblocks_tree,att_dim,error)
        call h5aclose_f(att_id,error)
        call h5aopen_f(file_id,'time',att_id,error)
        call h5aread_f(att_id,h5t_native_double,time_sys%t,att_dim,error)
        call h5aclose_f(att_id,error)
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(nblocks_tree,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(time_sys%t,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call read_tree(file_id)
    if (rank==0) then
        call h5fclose_f(file_id,error)
        call h5close_f(error)
        call chdir(path_root)
    end if
end subroutine hdf_read

subroutine read_tree(file_id)
    !read in the keys and rebuild the tree
    type(blockdef), pointer :: blk_traversal
    integer(HID_T) :: file_id,dset_id
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer :: i,error,ierr
    integer, dimension(:,:), allocatable :: keys
    allocate(dims(2),keys(nblocks_tree,6))
    if (rank==0) then
        dims=(/nblocks_tree,6/)
        call h5dopen_f(file_id,'block_tree',dset_id,error)
        call h5dread_f(dset_id,h5t_native_integer,keys,dims,error)
        call h5dclose_f(dset_id,error)
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(keys,nblocks_tree*6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call rebuild_tree(keys)
    call rebuild_linked_list()
    call smr_load_balancer()
    call identify_block_spatial_type()
    call reallocate_all_domain_blocks()
    call link_neighbour_blocks()
    call read_heavy_data(file_id)
    call distribute_heavy_data()
    call restart_complete_all_data()
    call applyboundconds()
    call restart_static_level()
    deallocate(dims,keys)
end subroutine read_tree

subroutine rebuild_tree(keys)
    !preorder traversal
    !logical level have level=-1, domain blocks have level>=0
    !rebuild the links within a node. No inter nodes links
    type(blockdef), pointer :: blk_traversal
    integer, dimension(:,:), allocatable :: keys
    integer :: i,j,key_p(3),key(3),level,it
    real(8) :: coords(2),n_size(2)
    call calculate_domain_block_parameters()
    if (nd==1) then
        call new_blank_block(blk_root)
        blk_root%key=keys(1,4:6)
        blk_root%nx_domain=nx
        blk_traversal=>blk_root
        do i=2,nblocks_tree
            key=keys(i,4:6)
            key_p=keys(i,1:3)
            if (key(3)>blk_traversal%key(3)) then
                call new_blank_block(blk_traversal%blk_xl)
                blk_traversal%blk_xl%blk_p=>blk_traversal
                blk_traversal%blk_xl%key=key
                blk_traversal=>blk_traversal%blk_xl
                if (key(3)>=llevel_max) then
                    level=key(3)-llevel_max
                    coords(1)=n_domain(1)+(n_domain(2)-n_domain(1))/nx_blks/2**level*(blk_traversal%key(1)-1)
                    coords(2)=0d0
                    n_size(1)=(n_domain(2)-n_domain(1))/nx_blks/2**level
                    n_size(2)=0d0
                    call new_block(blk_traversal,level,coords,n_size)
                else
                    blk_traversal%level=-1
                end if
            else
                it=1+blk_traversal%key(3)-key(3)
                do j=1,it
                    blk_traversal=>blk_traversal%blk_p
                end do
                call new_blank_block(blk_traversal%blk_xu)
                blk_traversal%blk_xl%blk_next=>blk_traversal%blk_xu
                blk_traversal%blk_xu%blk_pre=>blk_traversal%blk_xl
                blk_traversal%blk_xu%blk_p=>blk_traversal
                blk_traversal%blk_xu%key=key
                blk_traversal=>blk_traversal%blk_xu
                if (key(3)>=llevel_max) then
                    level=key(3)-llevel_max
                    coords(1)=n_domain(1)+(n_domain(2)-n_domain(1))/nx_blks/2**level*(blk_traversal%key(1)-1)
                    coords(2)=0d0
                    n_size(1)=(n_domain(2)-n_domain(1))/nx_blks/2**level
                    n_size(2)=0d0
                    call new_block(blk_traversal,level,coords,n_size)
                else
                    blk_traversal%level=-1
                end if
            end if
        end do
    else if (nd==2) then
    end if
end subroutine rebuild_tree

subroutine restart_static_level()
    !the statically refined region may not be the same as the one before restart
    !upon restart, the maximum level change of SMR is 1
    type(blockdef), pointer :: blk
    integer :: i,static_level,ierr,j
    real(8) :: zone(4)
    if (refine_type=='mixed'.or.refine_type=='static') then
        do i=1,nrefine_region
            call reset_smr_processor_variables()
            zone=refine_region(i,1:4)
            static_level=int(refine_region(i,5))
            call recalculate_static_level(blk_root,zone,static_level,processor%smr_ngrow)
            call mpi_gather(processor%smr_ngrow,1,MPI_INTEGER,processor%mpi_processor_oper,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(processor%mpi_processor_oper,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            np_nblk=np_nblk+processor%mpi_processor_oper
        end do
        call restart_static_buffer_blocks()
        !call restart_static_trim_blocks(blk_root)
        call restart_smr_load_balancer()
        call rank0_allocate()
        call identify_block_spatial_type()
        call unlink_neighbour_blocks(blk_root)
        call link_neighbour_blocks()
        call communicate_hydro()
    end if
end subroutine restart_static_level

subroutine recalculate_static_level(blk,zone,static_level,nrefine)
    !label the blocks that overlap the SMR zones, refine if the blk belongs is on the processor
    !maintain the linked list
    type(blockdef), pointer :: blk
    real(8) :: zone(4)
    integer :: static_level,nrefine
    if (nd==1) then
        if (associated(blk%blk_xl)) then
            !blk has children
            if (block_zone_overlap(blk,zone).or.blk%key(3)<llevel_max) then
                blk%static_level=blk%level
                call recalculate_static_level(blk%blk_xl,zone,static_level,nrefine)
                call recalculate_static_level(blk%blk_xu,zone,static_level,nrefine)
            end if
        else
            !blk does not have children
            if (block_zone_overlap(blk,zone).and.blk%level<static_level) then
                blk%static_level=blk%level
                call grow_static_tree_node(blk,zone,static_level)
                if (blk%on_processor) then
                    call recursive_allocate_heavy_data(blk,nrefine)
                end if
            end if
        end if
    else if (nd==2) then
    end if
end subroutine recalculate_static_level

subroutine recursive_allocate_heavy_data(blk,nrefine)
    type(blockdef), pointer :: blk
    integer :: nrefine
    if (associated(blk%blk_xl)) then
        call allocate_block_heavy_data(blk%blk_xl)
        call allocate_block_heavy_data(blk%blk_xu)
        call block_prolongation(blk)
        nrefine=nrefine+1
        call recursive_allocate_heavy_data(blk%blk_xl,nrefine)
        call recursive_allocate_heavy_data(blk%blk_xu,nrefine)
    end if
end subroutine recursive_allocate_heavy_data

subroutine restart_static_buffer_blocks()
    !check if every domain block's neighbour has buffer blocks in SMR
    integer :: i,level,ierr
    do i=1,max_refine_level-1
        call reset_smr_processor_variables()
        call restart_check_static_buffer_blocks(blk_root,processor%smr_ngrow)
        call mpi_gather(processor%smr_ngrow,1,MPI_INTEGER,processor%mpi_processor_oper,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(processor%mpi_processor_oper,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        np_nblk=np_nblk+processor%mpi_processor_oper
    end do
end subroutine restart_static_buffer_blocks

subroutine restart_check_static_buffer_blocks(blk,nrefine)
    type(blockdef), pointer :: blk,blk_temp
    integer :: key(3),i,level,nrefine
    if (nd==1) then
        if (associated(blk%blk_xl)) then
            !left neighbour
            if (blk%blk_xl%static_level>1.and.(.not.is_xl_block(blk%blk_xl))) then
                key=key_xl_neighbour(blk%blk_xl%key)
                key=key_parent(key)
                call block_key_to_pointer(key,blk_temp)
                if (blk%blk_xl%key(3)-blk_temp%key(3)>1) then
                    !neighbour's parent is not created
                    call grow_tree_node_1d(blk_temp)
                    if (blk_temp%on_processor) then
                        call allocate_block_heavy_data(blk_temp%blk_xl)
                        call allocate_block_heavy_data(blk_temp%blk_xu)
                        call block_prolongation(blk_temp)
                        print *,blk%w(blk_xlb:blk_xub,1,1,1),'haha'
                        print *,blk%blk_xl%w(blk_xlb:blk_xub,1,1,1)
                        print *,blk%blk_xu%w(blk_xlb:blk_xub,1,1,1)
                        stop
                        nrefine=nrefine+1
                    end if
                    level=blk_temp%level+1
                    blk_temp%static_level=blk_temp%level
                    blk_temp%blk_xl%static_level=level
                    blk_temp%blk_xu%static_level=level
                else
                    if (blk%blk_xl%static_level-blk_temp%static_level>1) then
                        blk_temp%static_level=blk%blk_xl%static_level-1
                    end if
                end if
            end if
            call restart_check_static_buffer_blocks(blk%blk_xl,nrefine)
        end if
        if (associated(blk%blk_xu)) then
            !right neighbour
            if (blk%blk_xu%static_level>1.and.(.not.is_xu_block(blk%blk_xu))) then
                key=key_xu_neighbour(blk%blk_xu%key)
                key=key_parent(key)
                call block_key_to_pointer(key,blk_temp)
                if (blk%blk_xu%key(3)-blk_temp%key(3)>1) then
                    !neighbour's parent is not created
                    call grow_tree_node_1d(blk_temp)
                    if (blk_temp%on_processor) then
                        call allocate_block_heavy_data(blk_temp%blk_xl)
                        call allocate_block_heavy_data(blk_temp%blk_xu)
                        call block_prolongation(blk_temp)
                        print *,blk%w(blk_xlb:blk_xub,1,1,1),'hoho'
                        print *,blk%blk_xl%w(blk_xlb:blk_xub,1,1,1)
                        print *,blk%blk_xu%w(blk_xlb:blk_xub,1,1,1)
                        nrefine=nrefine+1
                    end if
                    level=blk_temp%level+1
                    blk_temp%static_level=blk_temp%level
                    blk_temp%blk_xl%static_level=level
                    blk_temp%blk_xu%static_level=level
                else
                    if (blk%blk_xu%static_level-blk_temp%static_level>1) then
                        blk_temp%static_level=blk%blk_xu%static_level-1
                        if (associated(blk_temp%blk_next)) then
                            blk_temp%blk_next%static_level=blk%blk_xu%static_level-1
                        end if
                    end if
                end if
            end if
            call restart_check_static_buffer_blocks(blk%blk_xu,nrefine)
        end if
    else if (nd==2) then
    end if
end subroutine restart_check_static_buffer_blocks

subroutine restart_static_trim_blocks(blk)
    !if blk%static_level<blk%level derefine the blk
    !maintain the linked list
    !can only reduce 1 level at each restart
    type(blockdef), pointer :: blk,blk_xu,blk_xl
    integer :: i,level
    if (nd==1) then
        if (blk%level<=llevel_max) then
            !logic level, proceed
            if (associated(blk%blk_xl)) then
                blk_xl=>blk%blk_xl
                call restart_static_trim_blocks(blk_xl)
            end if
            if (associated(blk%blk_xu)) then
                blk_xu=>blk%blk_xu
                call restart_static_trim_blocks(blk_xu)
            end if
        else
            !refined level
            if (associated(blk%blk_xl)) then
                blk_xl=>blk%blk_xl
                blk_xu=>blk%blk_xu
                if (blk_xl%level>blk_xl%static_level.and.blk_xu%level>blk_xu%static_level) then
                    call block_restriction(blk)
                    if (associated(blk_xl,blk_processor_head)) then
                        blk_processor_head=>blk
                    end if
                    if (associated(blk_xu,blk_processor_tail)) then
                        blk_processor_tail=>blk
                    end if
                    if (associated(blk_xl%blk_pre)) then
                        blk_xl%blk_pre%blk_next=>blk
                        blk%blk_pre=>blk_xl%blk_pre
                    else            !this is the first block
                        blk_head=>blk
                        blk%blk_next=>blk_xu%blk_next
                        blk_xu%blk_next%blk_pre=>blk
                    end if
                    if (associated(blk_xu%blk_next)) then
                        blk%blk_next=>blk_xu%blk_next
                        blk_xu%blk_next%blk_pre=>blk
                    else            !this is the last block
                        blk_tail=>blk
                        blk%blk_pre=>blk_xl%blk_pre
                        blk_xl%blk_pre%blk_next=>blk
                    end if
                    call destroy_block(blk_xl)
                    call destroy_block(blk_xu)
                    nullify(blk%blk_xl,blk%blk_xu)
                else
                    call restart_static_trim_blocks(blk_xl)
                    call restart_static_trim_blocks(blk_xu)
                end if
            end if
        end if
    else if (nd==2) then
    end if
end subroutine restart_static_trim_blocks

subroutine rank0_allocate()
    type(blockdef), pointer :: blk
    integer :: i
    blk=>blk_head
    do i=1,nblk_total
        if (.not.allocated(blk%w)) then
            call allocate_block_heavy_data(blk)
        end if
        blk=>blk%blk_next
    end do
end subroutine rank0_allocate

subroutine rebuild_linked_list()
    !link the leaves only
    type(blockdef), pointer :: blk_builder1,blk_builder2,blk_builder3
    real(8) :: coords(2),n_size(2)
    integer :: i
    if (nd==1) then
        call find_the_first_node(blk_root,blk_head)
        blk_builder1=>blk_head
        do while (associated(blk_builder1))
            call find_the_successor_node(blk_builder1,blk_builder2)
            call find_the_first_node(blk_builder2,blk_builder3)
            if (.not.associated(blk_builder2,blk_root)) then
                blk_builder1%blk_next=>blk_builder3
                blk_builder3%blk_pre=>blk_builder1
                blk_builder1=>blk_builder3
            else
                exit
            end if
        end do
        blk_tail=>blk_builder1
        nullify(blk_builder1,blk_builder2)
    else if (nd==2) then
    end if
end subroutine rebuild_linked_list

subroutine reallocate_all_domain_blocks()
    type(blockdef), pointer :: blk_traversal
    blk_traversal=>blk_root
    call reallocate_blocks(blk_traversal)
end subroutine reallocate_all_domain_blocks

subroutine reallocate_blocks(blk)
    type(blockdef), pointer :: blk
    if (nd==1) then
        if (associated(blk%blk_xl)) then
            if (blk%blk_xl%key(3)>=llevel_max) then
                call allocate_block_heavy_data(blk%blk_xl)
            end if
            call reallocate_blocks(blk%blk_xl)
        end if
        if (associated(blk%blk_xu)) then
            if (blk%blk_xu%key(3)>=llevel_max) then
                call allocate_block_heavy_data(blk%blk_xu)
            end if
            call reallocate_blocks(blk%blk_xu)
        end if
    else if (nd==2) then
    end if
end subroutine reallocate_blocks

subroutine read_heavy_data(file_id)
    !master rank read all the data and store it to its blocks
    type(blockdef), pointer :: blk_traversal
    integer(HSIZE_T) :: dims(2)
    integer(HID_T) :: file_id,dset_id
    integer :: i,j,error,ndset,ierr
    character(len=16) :: fmt1,s1
    character(len=32) :: groupname,dsetname(100)
    if (rank==0) then
        blk_traversal=>blk_head
        dims=(/blk_size_nx,blk_size_ny/)
        fmt1='(I5.5)'
        if (nd==1) then
            ndset=4
            dsetname(1)='rho';dsetname(2)='vx';dsetname(3)='pres';dsetname(4)='temp'
            if (iradiation==4) then
                ndset=ndset+1
                dsetname(5)='Erad'
            end if
        else if (nd==2) then
        end if
        do i=1,nblk_total
            write(s1,fmt1) blk_traversal%blk_id
            groupname='blk'//trim(s1)
            do j=1,ndset
                call read_a_data(file_id,groupname,dsetname(j),blk_traversal,dims)
            end do
            blk_traversal=>blk_traversal%blk_next
        end do
        nullify(blk_traversal)
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
end subroutine read_heavy_data

subroutine distribute_heavy_data()
    !rank0 distribute the blocks to other blocks according to the load balancer
    type(blockdef), pointer :: blk
    integer :: blk_id,recv_rank,i,ierr
    if (np>1) then
        if (rank==0) then
            blk=>blk_processor_tail%blk_next
            do while (associated(blk))
                blk_id=blk%blk_id
                call block_id_to_processor_rank(np_nblk,blk_id,recv_rank)
                call send_block(blk,recv_rank)
                blk=>blk%blk_next
            end do
        end if
        if (rank/=0) then
            blk=>blk_processor_head
            do i=1,np_nblk(rank+1)
                blk_id=blk%blk_id
                call receive_block(blk,0)
                blk=>blk%blk_next
            end do
        end if
    end if
end subroutine distribute_heavy_data

subroutine read_a_data(file_id,groupname,dsetname,blk,dims)
    type(blockdef), pointer :: blk
    integer(HID_T) :: file_id,group_id,dset_id
    integer(HSIZE_T) :: dims(2)
    character(len=32) :: groupname,dsetname
    integer :: error
    call h5gopen_f(file_id,groupname,group_id,error)
    call h5dopen_f(group_id,trim(dsetname),dset_id,error)
    !call h5dopen_f(file_id,trim(blockname)//trim(dsetname),dset_id,error)
    select case (dsetname)
    case('rho')
        call h5dread_f(dset_id,h5t_native_double,blk%w(1:blk_size_nx,1:blk_size_ny,1,irho),dims,error)
    case('vx')
        call h5dread_f(dset_id,h5t_native_double,blk%w(1:blk_size_nx,1:blk_size_ny,1,ivx),dims,error)
    case('pres')
        call h5dread_f(dset_id,h5t_native_double,blk%w(1:blk_size_nx,1:blk_size_ny,1,ipres),dims,error)
    case('temp')
        call h5dread_f(dset_id,h5t_native_double,blk%temp(1:blk_size_nx,1:blk_size_ny,1),dims,error)
    case('Erad')
        call h5dread_f(dset_id,h5t_native_double,blk%Erad(1:blk_size_nx,1:blk_size_ny,1),dims,error)
    end select
    call h5dclose_f(dset_id,error)
    call h5gclose_f(group_id,error)
end subroutine read_a_data

subroutine restart_complete_all_data()
    !calculate the derived physical quantities
    !fill all guard cells of all blocks
    type(blockdef), pointer :: blk_traversal
    integer :: i,j,iblk
    real(8) :: rho,temp,egv,w(5),u(5)
    integer :: ierr
    blk_traversal=>blk_processor_head
    do j=1,np_nblk(rank+1)
        if (nd==1) then
            do i=1,blk_size_nx
                rho=blk_traversal%w(i,1,1,1)
                temp=blk_traversal%temp(i,1,1)
                w=blk_traversal%w(i,1,1,1:5)
#if             ieos==1
                egv=calculate_egv_from_rho_temp(rho,temp)
                call wtou(w,u)
#elif           ieos==2
                egv=egvrhot(rho,temp)
                call eos_hllc_analytic_wtou(w,egv,u)
#endif
                blk_traversal%egv(i,1,1)=egv
                blk_traversal%u(i,1,1,1:5)=u
            end do
        else if (nd==2) then
        end if
        blk_traversal=>blk_traversal%blk_next
    end do
    nullify(blk_traversal)
    call communicate_hydro()
end subroutine restart_complete_all_data

end module io_in
