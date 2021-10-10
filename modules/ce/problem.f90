module problem
use phylib
use datastructure
use mathlib
use eos
implicit none
real(8), dimension(5), protected :: w0,u0
real(8), protected :: rho0,v0(3),temp0,p0,egv0,E0,dx,m_star,v_esc,  &
    t_outflow1,temp_outflow1,mdot1,v_esc_coefficient1,p_outflow1,rho_outflow1,v_outflow1,egv_outflow1,ek_outflow1,hii_outflow1, &
    t_outflow2,temp_outflow2,mdot2,v_esc_coefficient2,p_outflow2,rho_outflow2,v_outflow2,egv_outflow2,ek_outflow2,hii_outflow2, &
    dt_record,ek_rank,eg_rank,rho_crit,t_interval,   &
    eg,ek,eg_pre,ek_pre,ek_total,m_total,rtau,rho_floor,        &
    eg_total,degdt,dekdt,eg_flux,e_flux,gp_rank,gp,gp_pre,dgpdt,dt_pre,er,er_rank,er_pre,derdt
logical :: loc_tau
logical, dimension(:), allocatable :: mpi_tau


contains

subroutine initialize_problem()
    real(8), allocatable :: species(:)
    !calculate w, u, temp, and egv for all location
    namelist /rhd_quantities/ out_dir,rho0,v0,temp0,m_star,t_outflow1,temp_outflow1,mdot1,v_esc_coefficient1,    &
    t_outflow2,temp_outflow2,mdot2,v_esc_coefficient2,dt_record,petsc_rtol,record_length,energy_conserve_formalism,     &
    petsc_eos_rtol,resolve_rad_dt,erad_slope_refine_thresh,rtau,rho_crit,t_interval,opacity_gas_rho_min,opacity_gas_rho_max,   &
    opacity_gas_t_min,opacity_gas_t_max,rho_floor
    open(unit=11,file='modules/problem/problem.data',status='old',action='read')
    read(unit=11,nml=rhd_quantities)
    path_out=trim(path_root)//'/modules/problem/'//trim(out_dir)
    close(11)
    allocate(mpi_tau(np))
    eg_pre=0d0
    ek_pre=0d0
    gp_pre=0d0
    t_outflow1=t_outflow1*day
    t_outflow2=t_outflow2*day
    t_interval=t_interval*day
    w0(1)=rho0
    w0(2:4)=v0
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
    E0=a_rad*temp0**4
    dx=dxyz(1)
    central_star%core%mass=m_star
    v_esc=sqrt(2*gr*m_star/n_domain(1))
    mdot1=mdot1*msun/year
    mdot2=mdot2*msun/year
    time_sys%dt_record=dt_record
    v_outflow1=v_esc*v_esc_coefficient1
    v_outflow2=v_esc*v_esc_coefficient2
    rho_outflow1=mdot1/4/pi/n_domain(1)**2/v_outflow1
    rho_outflow2=mdot1/4/pi/n_domain(1)**2/v_outflow2
#if     ieos==1
    p_outflow1=idealgas_prhot(rho_outflow1,temp_outflow1)
    p_outflow2=idealgas_prhot(rho_outflow2,temp_outflow2)
    egv_outflow1=idealgas_egv(p_outflow1)
    egv_outflow2=idealgas_egv(p_outflow2)
#elif   ieos==2
    p_outflow1=prhot(rho_outflow1,temp_outflow1)
    p_outflow2=prhot(rho_outflow2,temp_outflow2)
    egv_outflow1=egvrhot(rho_outflow1,temp_outflow1)
    egv_outflow2=egvrhot(rho_outflow2,temp_outflow2)
    allocate(species(5))
    call calspecies(rho_outflow1,temp_outflow1,species)
    hii_outflow1=species(3)/(2*species(1)+species(2)+species(3))
    call calspecies(rho_outflow2,temp_outflow2,species)
    hii_outflow2=species(3)/(2*species(1)+species(2)+species(3))
    deallocate(species)
#endif
    ek_total=mdot1*t_outflow1*half*v_outflow1*v_outflow1+mdot2*t_outflow2*half*v_outflow2*v_outflow2
    eg_total=mdot1*t_outflow1*egv_outflow1/rho_outflow1+mdot2*t_outflow2*egv_outflow2/rho_outflow2
    m_total=(mdot1*t_outflow1+mdot2*t_outflow2)/msun
    !initialize all the cells including boundary cells
    if (rank==0) then
        write(*,'(A32,7ES18.6E2)') 'primitive quantities',w0
        write(*,'(A32,7ES18.6E2)') 'conserved quantities',u0
        write(*,'(A32,ES18.6E2)') 'escape velocity',v_esc
        write(*,'(A32,ES18.6E2)') 'dt_record',time_sys%dt_record
        call chdir(path_out)
        open(unit=15,file='setup.dat',status='replace',action='write')
        write(unit=15,fmt=*) '&setup'
        write(unit=15,fmt=*) '  m_total=',m_total
        write(unit=15,fmt=*) '  ek_total=',ek_total
        write(unit=15,fmt=*) '  eg_total=',eg_total
        write(unit=15,fmt=*) '  m1=',mdot1*t_outflow1/msun
        write(unit=15,fmt=*) '  rho_outflow1=',rho_outflow1
        write(unit=15,fmt=*) '  v_outflow1=',v_outflow1
        write(unit=15,fmt=*) '  temp_outflow1=',temp_outflow1
        write(unit=15,fmt=*) '  hii_fraction1=',hii_outflow1
        write(unit=15,fmt=*) '  m2=',mdot2*t_outflow2/msun
        write(unit=15,fmt=*) '  rho_outflow2=',rho_outflow2
        write(unit=15,fmt=*) '  v_outflow2=',v_outflow2
        write(unit=15,fmt=*) '  temp_outflow2=',temp_outflow2
        write(unit=15,fmt=*) '  hii_fraction2=',hii_outflow2
        write(unit=15,fmt=*) '/'
        close(15)
        call chdir(path_root)
    end if
end subroutine initialize_problem

subroutine initial_hydro(pos,t,w,u,temp,egv,mark)
    real(8) :: pos(3),t,w(5),u(5),temp,egv,r,r0,rr
    integer, optional :: mark
    r=norm2(pos)
    r0=n_domain(1)
    rr=r/r0
    w(1)=max(w0(1)/rr**2*sqrt(rr),rho_floor)
    !w(1)=max(w0(1)/rr**2,rho_floor)
    w(2)=w0(2)/sqrt(rr)
    w(3:4)=0d0
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
    dw=zero
end subroutine perturb

subroutine time_dependent_bound_type()
    real(8) :: t
    t=time_sys%t
    if (t>t_outflow1.and.t<(t_outflow1+t_interval)) then
        hydro_bound_type(1)=1
    else if (t>(t_outflow1+t_outflow2+t_interval)) then
        hydro_bound_type(1)=1
    else
        hydro_bound_type(1)=9
    end if
end subroutine time_dependent_bound_type

subroutine boundary_hydro(pos,t,w,u,temp,egv,mark)
    real(8) :: pos(3),t,w(5),u(5),temp,egv,p,rho,v
    real(8) :: mdot,v_outflow
    integer, optional :: mark
    if (present(mark)) then
        if (mark==1) then
            if (t<=t_outflow1) then
                w(1)=rho_outflow1
                w(2)=v_outflow1
                w(3:4)=0d0
                temp=temp_outflow1
#if         ieos==1
                w(5)=p_outflow1
                call wtou(w,u)
#elif       ieos==2
                !specify the s of the second ejecta
                egv=egvrhot(rho_outflow1,temp_outflow1)
                w(5)=prhot(rho_outflow1,temp_outflow1)
                call eos_hllc_analytic_wtou(w,egv,u)
#endif
            else if (t>t_outflow1+t_interval.and.t<=(t_outflow1+t_outflow2+t_interval)) then
                w(1)=rho_outflow2
                w(2)=v_outflow2
                w(3:4)=0d0
                temp=temp_outflow2
#if         ieos==1
                w(5)=p_outflow2
                call wtou(w,u)
#elif       ieos==2
                !specify the s of the second ejecta
                egv=egvrhot(rho_outflow2,temp_outflow2)
                w(5)=prhot(rho_outflow2,temp_outflow2)
                call eos_hllc_analytic_wtou(w,egv,u)
#endif
            end if
!            if (t<=t_outflow1+t_outflow2) then
!                mdot=mdot1+(mdot2-mdot1)*t/(t_outflow1+t_outflow2)
!                v_outflow=v_outflow1+(v_outflow2-v_outflow1)*t/(t_outflow1+t_outflow2)
!                temp=temp_outflow1+(temp_outflow2-temp_outflow1)*t/(t_outflow1+t_outflow2)
!                rho=mdot/4/pi/n_domain(1)**2/v_outflow
!                egv=egvrhot(rho,temp)
!                p=prhot(rho,temp)
!                w(1)=rho
!                w(2)=v_outflow
!                w(3:4)=0d0
!                w(5)=p
!                call eos_hllc_analytic_wtou(w,egv,u)
!            end if
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
    if (rank==np-1) then
        blk=>blk_tail
        if (iradiation==4) then
            record_array(1)=4d0*pi*n_domain(2)**2*blk%Fradx(blk_size_nx,1,1)
        else
            record_array(1)=0d0!4d0*pi*n_domain(2)**2*blk%Fradx(blk_size_nx,1,1)/lsun
        end if
    end if
    call mpi_bcast(record_array,1,MPI_REAL8,np-1,MPI_COMM_WORLD,ierr)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    record_array(2)=eg/lsun/day
    record_array(3)=ek/lsun/day
    record_array(4)=er/lsun/day
    record_array(5)=degdt/lsun
    record_array(6)=dekdt/lsun
    record_array(7)=derdt/lsun
    if (rank==0) then
        blk=>blk_head
        r=n_domain(1)
        eg_flux=blk%egv(0,1,1)*4d0*pi*r**2*blk%w(0,1,1,2)
        if (blk%xflux(0,1,1,1)>0) then
            phi=blk%gpotential(0,1,1)
        else
            phi=blk%gpotential(1,1,1)
        end if
        e_flux=(blk%xflux(0,1,1,5)+phi*blk%xflux(0,1,1,1))*4d0*pi*r**2
        record_array(9)=e_flux/lsun
        !record_array(10)=blk%xflux(0,1,1,1)!blk%w(1,1,1,1)*4d0*pi*r**2*blk%w(1,1,1,2)*year/msun
    end if
    record_array(8)=dgpdt/lsun
    record_array(10)=record_array(5)+record_array(6)+record_array(7)+record_array(8)-record_array(9)
    !record_array(11)=record_array(5)+record_array(6)+record_array(8)-record_array(9)
    record_array(11)=time_sys%t
    !r_tau=rosseland_tau_radius(rtau)
    r_tau=ejecta_radius(rho_crit)
    record_array(12)=pow(record_array(1)/sigma_sb/4d0/pi/r_tau**2,0.25d0)
    record_array(13)=r_tau/rsun
end subroutine assemble_record_array

subroutine problem_oper()
    integer :: ierr
    if (time_sys%ntimestep==0) then
        if (rank==0) write(*,fmt='(32A14)') 'luminosity','eg','ek','er','degdt','dekdt','derdt','dgpdt',    &
            '(E+p+phi)flux','detotaldt','t','Teff','R'
    end if
    eg_rank=sum_eg()
    ek_rank=sum_ek()
    er_rank=sum_er()
    gp_rank=sum_gpotential()
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_allreduce(eg_rank,eg,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_allreduce(ek_rank,ek,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_allreduce(er_rank,er,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_allreduce(gp_rank,gp,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    degdt=(eg-eg_pre)/time_sys%dt_hydro
    dekdt=(ek-ek_pre)/time_sys%dt_hydro
    derdt=(er-er_pre)/time_sys%dt_hydro
    dgpdt=(gp-gp_pre)/time_sys%dt_hydro
    eg_pre=eg
    ek_pre=ek
    er_pre=er
    gp_pre=gp
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

function ejecta_radius(rho_crit)
    type(blockdef), pointer :: blk
    real(8) :: ejecta_radius,rho_crit,r
    integer :: i,j,ierr
    processor%processor_logical=.false.
    processor%global_logical=.false.
    blk=>blk_processor_tail
    outer: do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            if (blk%w(blk_size_nx+1-j,1,1,1)>rho_crit) then
                processor%processor_logical=.true.
                r=blk%x_center(blk_size_nx+1-j)
                exit outer
            end if
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_pre
        end if
    end do outer
    call mpi_gather(processor%processor_logical,1,MPI_LOGICAL,processor%global_logical,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(processor%global_logical,np,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    do i=np,1,-1
        if (processor%global_logical(i)) then
            call mpi_bcast(r,1,MPI_REAL8,i-1,MPI_COMM_WORLD,ierr)
            exit
        end if
    end do
    ejecta_radius=r
end function ejecta_radius

subroutine finalize_problem()
end subroutine finalize_problem

end module problem
