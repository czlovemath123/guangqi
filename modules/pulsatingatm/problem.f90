module problem
use phylib
use datastructure
use mathlib
use eos
implicit none
real(8), dimension(5), protected :: w0,u0
real(8), protected :: m_star,rho0,v0,temp0,dt_record,period,egv0,p0,v_esc
logical :: loc_tau
logical, dimension(:), allocatable :: mpi_tau


contains

subroutine initialize_problem()
    real(8), allocatable :: species(:)
    !calculate w, u, temp, and egv for all location
    namelist /rhd_quantities/ m_star,period,dt_record,rho0,v0,temp0,petsc_rtol
    open(unit=11,file='modules/problem/problem.data',status='old',action='read')
    read(unit=11,nml=rhd_quantities)
    path_out=trim(path_root)//'/modules/problem/'//trim(out_dir)
    close(11)
    allocate(mpi_tau(np))
    w0(1)=rho0
    w0(2:4)=0d0
#if     ieos==1
    p0=idealgas_prhot(rho0,temp0)
    w0(5)=p0
    egv0=idealgas_egv(p0)
    call wtou(w0,u0)
#elif   ieos==2
    p0=prhot(rho0,temp0)
    w0(5)=p0
    egv0=egvrhot(rho0,temp0)
    call eos_hllc_analytic_wtou(w0,egv0,u0)
#endif
    central_star%core%mass=m_star
    v_esc=sqrt(2*gr*m_star/n_domain(1))
    time_sys%dt_record=dt_record
    !initialize all the cells including boundary cells
    if (rank==0) then
        write(*,'(A32,7ES18.6E2)') 'primitive quantities',w0
        write(*,'(A32,7ES18.6E2)') 'conserved quantities',u0
        write(*,'(A32,ES18.6E2)') 'escape velocity',v_esc
        write(*,'(A32,ES18.6E2)') 'dt_record',time_sys%dt_record
        call chdir(path_out)
        open(unit=15,file='setup.dat',status='replace',action='write')
        close(15)
        call chdir(path_root)
    end if
end subroutine initialize_problem

subroutine initial_hydro(pos,t,w,u,temp,egv,mark)
    real(8) :: pos(3),t,w(5),u(5),temp,egv,r,r0,rr
    integer, optional :: mark
    r=norm2(pos)
    r0=n_domain(1)
    rr=40*r/r0
    w(1)=w0(1)/rr**2/sqrt(rr)
    w(2:4)=0d0
    temp=temp0/rr
#if     ieos==1
    w(5)=idealgas_prhot(w(1),temp)
    egv=idealgas_egv(w(5))
    call wtou(w,u)
#elif   ieos==2
    w(5)=prhot(w(1),temp)
    egv=egvrhot(w(1),temp)
    call eos_hllc_analytic_wtou(w,egv,u)
#endif
end subroutine initial_hydro

subroutine initial_rad(pos,t,Erad,Frad,feddington,mark)
    real(8) :: pos(3),t,Frad(3),feddington,r,r0,rr,temp
    real(8), dimension(:), allocatable :: Erad
    integer, optional :: mark
    r=norm2(pos)
    r0=n_domain(1)
    rr=r/r0
    temp=temp0/rr
    Erad=a_rad*temp**4
end subroutine initial_rad

subroutine perturb(pos,t,w,dw)
    real(8) :: pos(3),t,w(5),dw(5)
end subroutine perturb

subroutine time_dependent_bound_type()
end subroutine time_dependent_bound_type

subroutine boundary_hydro(pos,t,w,u,temp,egv,mark)
    real(8) :: pos(3),t,w(5),u(5),temp,egv,p,rho,v
    real(8) :: mdot,v_outflow
    integer, optional :: mark
    if (present(mark)) then
        if (mark==1) then
            w(1)=rho0
            w(2)=v0*sin(2d0*pi*t/period)
            w(3:4)=0d0
            temp=temp0
#if         ieos==1
            w(5)=p0
            call wtou(w,u)
#elif       ieos==2
            !specify the s of the second ejecta
            egv=egvrhot(rho0,temp0)
            w(5)=prhot(rho0,temp0)
            call eos_hllc_analytic_wtou(w,egv,u)
#endif
        else if (mark==2) then
            w(1)=1d-18!blk_tail%w(blk_size_nx,1,1,1)
            w(2)=0d0
            w(3:4)=0d0
#if         ieos==1
            temp=blk_tail%temp(blk_size_nx,1,1)
            w(5)=idealgas_prhot(w(1),temp)
            call wtou(w,u)
#elif       ieos==2
            !specify the s of the second ejecta
            egv=egvrhot(rho0,temp0)
            w(5)=prhot(rho0,temp0)
            call eos_hllc_analytic_wtou(w,egv,u)
#endif
        end if
    end if
end subroutine boundary_hydro

subroutine boundary_rad(pos,t,Erad,Frad,feddington,mark)
    real(8) :: pos(3),t,Frad(3),feddington,temp
    real(8), dimension(:), allocatable :: Erad
    integer, optional :: mark
end subroutine boundary_rad

subroutine assemble_record_array(record_array)
    type(blockdef), pointer :: blk
    real(8), allocatable :: record_array(:)
    real(8) :: r,m,phi,r_tau
    integer :: i,ierr
end subroutine assemble_record_array

subroutine problem_oper()
    integer :: ierr
end subroutine problem_oper

function rosseland_tau_radius(tau)
    !find the radius of tau
    type(blockdef), pointer :: blk
    real(8) :: rosseland_tau_radius,tau_rank_end,tau_rank_start,tau
    integer :: i,j,ierr,stat(MPI_STATUS_SIZE),loc_rank
    tau_rank_start=0d0
    blk=>blk_processor_tail
    do i=1,np_nblk(rank+1)
        if (associated(blk,blk_processor_tail)) then
            blk%tau(blk_size_nx,1,1)=0d0
        else
            blk%tau(blk_size_nx,1,1)=blk%blk_next%tau(0,1,1)
        end if
        do j=1,blk_size_nx
            blk%tau(blk_size_nx-j,1,1)=blk%tau(blk_size_nx-j+1,1,1)+blk%sigma_planck(blk_size_nx-j+1,1,1)*blk%dxyz(1)
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_pre
        end if
    end do
    tau_rank_start=blk_processor_tail%tau(blk_size_nx,1,1)
    tau_rank_end=blk_processor_head%tau(0,1,1)
    if (np>1) then
        if (rank==np-1) then
            call mpi_send(tau_rank_end,1,MPI_REAL8,rank-1,1,MPI_COMM_WORLD,ierr)
        else if (rank==0) then
            call mpi_recv(tau_rank_start,1,MPI_REAL8,rank+1,1,MPI_COMM_WORLD,stat,ierr)
            blk=>blk_processor_tail
            do i=1,np_nblk(rank+1)
                blk%tau=blk%tau+tau_rank_start
                if (i/=np_nblk(rank+1)) then
                    blk=>blk%blk_pre
                end if
            end do
        else
            call mpi_recv(tau_rank_start,1,MPI_REAL8,rank+1,1,MPI_COMM_WORLD,stat,ierr)
            blk=>blk_processor_tail
            do i=1,np_nblk(rank+1)
                blk%tau=blk%tau+tau_rank_start
                if (i/=np_nblk(rank+1)) then
                    blk=>blk%blk_pre
                end if
            end do
            call mpi_send(tau_rank_end,1,MPI_REAL8,rank-1,1,MPI_COMM_WORLD,ierr)
        end if
        tau_rank_start=blk_processor_tail%tau(blk_size_nx,1,1)
        tau_rank_end=blk_processor_head%tau(0,1,1)
    end if
    mpi_tau=.false.;loc_tau=.false.
    if (tau>tau_rank_start.and.tau<tau_rank_end) then
        loc_tau=.true.
    end if
    call mpi_gather(loc_tau,1,MPI_LOGICAL,mpi_tau,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(mpi_tau,np,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    !if (rank==0) print *,mpi_tau,tau_rank_start,tau_rank_end
    if (any(mpi_tau)) then
        do j=1,np
            if (mpi_tau(j)) then
                loc_rank=j-1
                exit
            end if
        end do
        if (loc_tau) then
            blk=>blk_processor_tail
            outer: do i=1,np_nblk(rank+1)
                do j=0,blk_size_nx
                    if (blk%tau(blk_size_nx-j,1,1)>=tau) then
                        rosseland_tau_radius=blk%pos(1)+blk%dxyz(1)*(blk_size_nx-j)
                        exit outer
                    end if
                end do
                if (i/=np_nblk(rank+1)) then
                    blk=>blk%blk_pre
                end if
            end do outer
        end if
        call mpi_bcast(rosseland_tau_radius,1,MPI_REAL8,loc_rank,MPI_COMM_WORLD,ierr)
    else
        rosseland_tau_radius=n_domain(1)
    end if
end function rosseland_tau_radius

subroutine finalize_problem()
end subroutine finalize_problem

end module problem
