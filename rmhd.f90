program rmhd
use hdf5
use mathlib
use datastructure
use eos
use communication
use phylib
use boundary
use io_in
use io_out
use gravity
use hydroscheme
use problem
use radiation
use cooling
use geometricsource
use source_control
use vis
implicit none

integer :: c1,c2,count_rate,count_max
integer :: ierr,errcode
real(8) :: dt_phy,t1,t2,dt_cpu
character(len=32) :: settings(3)
call problem_gen()
call problem_start()
call cpu_time(t1)
call system_clock(c1,count_rate,count_max)
call advan()
call system_clock(c2,count_rate,count_max)
call cpu_time(t2)
dt_phy=real(c2-c1,kind=8)/count_rate
dt_cpu=t2-t1
if (rank==0) then
    print *,'physical time is:',dt_phy
    print *,'cpu time is:',dt_cpu
end if
call finalize()

contains

subroutine problem_gen()
    integer :: i,ierr,errcode
    type(blockdef), pointer :: blk
    nullify(initialcond_scalar,boundcond_scalar)
    call start_mpi()
    call initialize_simulation_parameters()
    call initialize_eos_environment()
    call calculate_llevel_max()
    call initialize_problem()
    call dsetname_filename()
    call simulation_config()
    if (restart) then
        call hdf_read()
        call block_environment_setup()
        time_sys%t_output_initial=time_sys%t
        time_sys%dt_frame=(time_sys%t_final-time_sys%t_output_initial)/(nframe-iframe)
        allocate(t_target_output(nframe-iframe))
        do i=1,nframe-iframe
            t_target_output(i)=time_sys%t_output_initial+time_sys%dt_frame*i
        end do
        time_sys%t_next=max(time_sys%t_output_initial+time_sys%dt_frame,time_sys%t_final)
        !call display_blocks()
        !stop
    else
        call build_logical_tree()
        call build_linked_base_blocks()
        call refine_static_tree_and_blocks()
        call smr_load_balancer()
        call identify_block_spatial_type()
        call allocate_all_block_heavy_data()
        call block_environment_setup()
        call link_neighbour_blocks()
        call assemble_communication_pattern()
        time_sys%t=zero
        time_sys%t_output_initial=zero
        time_sys%dt_frame=(time_sys%t_final-time_sys%t_output_initial)/nframe
        allocate(t_target_output(nframe))
        do i=1,nframe
            t_target_output(i)=time_sys%dt_frame*i
        end do
        time_sys%t_next=max(time_sys%t_output_initial+time_sys%dt_frame,time_sys%t_final)
        !call display_blocks()
        !stop
    end if
    call initialize_radiation_environment()
end subroutine problem_gen

subroutine problem_start()
    if (.not.restart) then
        call apply_hydro_condition(initial_hydro)
        call apply_rad_condition(initial_rad)
        if (igeometry==2) then
            call apply_scalar_condition(initialcond_scalar,spherical2d_vphi)
        end if
        call output_blocks()
    end if
end subroutine problem_start

subroutine first_step()
    if (nd==1) then
        call grow_amr_tree()
    else if (nd==2) then
    else
    end if
end subroutine first_step

!evolve with time. Hydro, radiation, magnetic field, geometry, and gravity are doing their job here.
subroutine advan()
    integer :: i,j,k,ierr
    real(8) :: u(5),w(5),temp,egv
    real(8), allocatable :: record_array(:)
    character(len=30) :: fmt1
    character(len=5) :: mark
    character(len=128) :: alert
    i=0;j=0
    call first_step()
    do while(time_sys%t<time_sys%t_final)
        iframe=iframe+1
        i=i+1
        !time_sys%t_next=min(time_sys%t+time_sys%dt_frame,time_sys%t_final)
        time_sys%t_next=t_target_output(i)
        if (rank==0) print *,'frame: ',iframe,'next output time: ',time_sys%t_next/timescale,'t_final: ',time_sys%t_final/timescale
        do while(time_sys%t<time_sys%t_next)
            call grow_amr_tree()
            call hydro_step()
            if (source) then
                call communicate_hydro()
                call applyboundconds()
                call blk_traversal(source_apply)
            end if
            if (viscous) then
                call communicate_hydro()
                call applyboundconds()
                call blk_traversal(viscous_hydro)
            end if
            call radiation_transfer()
            call problem_oper()
            time_sys%t=time_sys%t+time_sys%dt_hydro
#if     irecord==1
            if (time_sys%dt_record==0) then
                alert='dt_record=0'
                call abort_achilles(alert)
            else
                if (time_sys%t/time_sys%dt_record>j_record.and.time_sys%ntimestep>0) then
                    j_record=j_record+1
                    allocate(record_array(record_length))
                    call assemble_record_array(record_array)
                    call record(record_array,j_record)
                    deallocate(record_array)
                end if
            end if
#endif
            !if (rank==0) print *,time_sys%t,time_sys%dt_hydro,blk_head%w(1,1,1,1:5)
            time_sys%ntimestep=time_sys%ntimestep+1
        end do
        call output_blocks()
    end do
#if     irecord==1
    call save_recorded_data(j_record)
#endif
end subroutine advan

subroutine finalize()
    !deallocate all the arrays and do some final output
    integer :: ierr
    call finalize_problem()
    call finalize_radiation()
    !call finalize_hydro(info)
    !call finalize_source(info)
    call finalize_mpi()
end subroutine finalize

end program rmhd
