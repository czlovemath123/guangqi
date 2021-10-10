module petsc_fld
#include <petsc/finclude/petscksp.h>
use datastructure
use eos
use communication
use phylib
use boundary
use radiation_common_functions
use petsc_fld_2d
use petscksp
implicit none

real(8), dimension(:), allocatable, protected :: petsc_local_T1,petsc_local_T2,petsc_local_E1,petsc_local_E2,   &
    petsc_diffT,petsc_diffE,beta_rank
logical, protected :: petsc_not_converged
integer, protected :: n_local_T,n_local_E

contains 

subroutine petsc_fld_initialize(userctx,ierr)
    PetscErrorCode ierr
    type(axbsystem) userctx
    Mat     A
    Vec     b,x
    KSP    ksp
    PetscInt Ntot,submat_nrow,n_local_rows
    PetscInt, allocatable :: nc(:),d_nz(:),d_nnz(:),o_nz(:),o_nnz(:)
    integer :: i
    if (nd==1) then
        Ntot=2*nblk_total*blk_size_nx
        n_local_rows=2*np_nblk(rank+1)*blk_size_nx
        n_local_T=np_nblk(rank+1)*blk_size_nx
        n_local_E=np_nblk(rank+1)*blk_size_nx
        call MatCreateAIJ(PETSC_COMM_WORLD,n_local_rows,n_local_rows,Ntot,Ntot,4,PETSC_NULL_INTEGER,1,PETSC_NULL_INTEGER,A,ierr);CHKERRQ(ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,n_local_rows,Ntot,b,ierr);CHKERRQ(ierr)
        call VecDuplicate(b,x,ierr);CHKERRQ(ierr)
        call KSPCreate(PETSC_COMM_WORLD,ksp,ierr);CHKERRQ(ierr)
        call KSPSetOperators(ksp,A,A,ierr);CHKERRQ(ierr)
#if     ieos==2
        allocate(petsc_local_T1(n_local_T),petsc_local_E1(n_local_E),petsc_local_T2(n_local_T),petsc_local_E2(n_local_E))
        allocate(petsc_diffT(n_local_T),petsc_diffE(n_local_E),beta_rank(n_local_T))
#endif
        userctx%x = x
        userctx%b = b
        userctx%A = A
        userctx%ksp = ksp
    else if (nd==2) then
        Ntot=2*nblk_total*blk_size_nx*blk_size_ny
        n_local_rows=2*np_nblk(rank+1)*blk_size_nx*blk_size_ny
        call MatCreateAIJ(PETSC_COMM_WORLD,n_local_rows,n_local_rows,Ntot,Ntot,6,PETSC_NULL_INTEGER,2,PETSC_NULL_INTEGER,A,ierr);CHKERRQ(ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,n_local_rows,Ntot,b,ierr);CHKERRQ(ierr)
        call VecDuplicate(b,x,ierr);CHKERRQ(ierr)
        call KSPCreate(PETSC_COMM_WORLD,ksp,ierr);CHKERRQ(ierr)
        call KSPSetOperators(ksp,A,A,ierr);CHKERRQ(ierr)
        userctx%x = x
        userctx%b = b
        userctx%A = A
        userctx%ksp = ksp
    end if
end subroutine petsc_fld_initialize

subroutine petsc_fld_estimate_dt(dtemp)
    type(blockdef), pointer :: blk
    real(8) :: dtemp,dt
    real(8), dimension(:), allocatable :: dt_rad_block
    integer :: i,ierr
    logical :: cv_temp
    cv_temp=.true.
    time_sys%dt_radiation=0d0
    allocate(dt_rad_block(np_nblk(rank+1)))
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        call petsc_fld_dt_block(blk,dtemp,dt)
        dt_rad_block(i)=dt
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    dt=minval(dt_rad_block)
    deallocate(dt_rad_block)
    time_sys%dt_processor=dt
    call mpi_gather(time_sys%dt_processor,1,MPI_REAL8,dt_processors,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    if (rank==0) time_sys%dt_radiation=minval(dt_processors)
    call mpi_bcast(time_sys%dt_radiation,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end subroutine petsc_fld_estimate_dt

subroutine petsc_fld_dt_block(blk,dtemp,dt)
    type(blockdef), pointer :: blk
    real(8) :: dtemp,dt,cv,kappa_planck,rho,Erad,temp
    real(8), dimension(:), allocatable :: dt_rad
    integer :: i
    allocate(dt_rad(blk_size_nx))
    do i=1,blk_size_nx
        cv=blk%cv(i,1,1)
        kappa_planck=blk%kappa_planck(i,1,1)
        Erad=blk%Erad(i,1,1)
        rho=blk%w(i,1,1,1)
        temp=blk%temp(i,1,1)
        if (abs(a_rad*temp**4-Erad)/Erad<1e-5) then
            dt_rad(i)=t_hubble
        else
            !dt_rad(i)=cv*dtemp/kappa_planck/rho/c_light/abs(a_rad*temp**4-Erad)
            call petsc_fld_dt_cell(rho,temp,Erad,kappa_planck,cv,dt_rad(i))
        end if
    end do
    dt=minval(dt_rad)
    deallocate(dt_rad)
end subroutine petsc_fld_dt_block

subroutine petsc_fld_dt_cell(rho,t,erad,kappa_planck,cv,dt)
    !for realistic eos
    real(8) :: rho,t,erad,kappa_planck,dt,t1,t2,t3,t4,eg_power,powerdiff,cv,g,dtemp
    dtemp=50d0
    eg_power=a_rad*t**4
    powerdiff=eg_power-erad
    g=kappa_planck*rho*c_light*abs(powerdiff)
#if     ieos==2
#if     ieosmodule==1
    !t1, t2, t3, and t4 are the temperatures that separate h2 disassociate start, h2 disassociation complete,
    !ionization start, and ionization complete
    call divide_temperature_domain_known(rho,t1,t2,t3,t4)
    if (t<=t1) then
        if (powerdiff>0) then
            dt=t_hubble
        else
            dt=cv*max((t1-t),dtemp)/g
        end if
    else if (t>=t4) then
        if (powerdiff>0) then
            dt=cv*max((t-t4),petsc_eos_rtol*t)/g
        else
            dt=t_hubble
        end if
    else
        !dt=cv*t*1d-1/g
        if (t2<t3) then
            if (t<t3) then
                dt=cv*max(min(t-t2,t3-t),petsc_eos_rtol*t)/g
            else
                dt=cv*max(min(t-t2,t3-t),petsc_eos_rtol*t)/g
            end if
        else
            if (t<t2) then
                dt=cv*t*petsc_eos_rtol/g
            else
                dt=cv*t*petsc_eos_rtol/g
            end if
        end if
    end if
#elif   ieosmodule==2
    call divide_temperature_domain_known(rho,t1,t2)
    if (t<t1) then
        if (powerdiff>0) then   !cooling
            dt=t_hubble
        else                    !heating
            dt=cv*max((t1-t),dtemp)/kappa_planck/rho/c_light/abs(powerdiff)
        end if
    else if (t>=t1.and.t<=t2) then
        dt=cv*max((t1-t),dtemp)/kappa_planck/rho/c_light/abs(powerdiff)
    else
        if (powerdiff>0) then   !cooling
            dt=cv*max((t-t2),dtemp)/kappa_planck/rho/c_light/abs(powerdiff)
        else                    !heating
            dt=t_hubble
        end if
    end if
#endif
#else
    dt=cv*t*1d-1/g
#endif
end subroutine petsc_fld_dt_cell

subroutine petsc_fld_solve(userctx,ierr)
    type(blockdef), pointer :: blk
    PetscErrorCode ierr
    type(axbsystem) userctx
    PetscScalar, pointer :: g(:)
    PetscInt, allocatable, dimension(:) :: mat_blks
    KSP, allocatable, dimension(:) :: subksp
    real(8) :: dt,dx,v,petsc_rtol_temp,dtemp,ddt
    PC   pc,subpc
    KSP  ksp
    Vec  b,x
    Mat  A
    integer :: Ntot,i,j,k,II,JJ,nlocal,first,nblocks,nsubcycle,blk_id,m,n,iter
    x    = userctx%x
    b    = userctx%b
    A    = userctx%A
    ksp  = userctx%ksp
    call communicate_hydro()
    call applyboundconds()
    if (nd==1) then
        call blk_traversal(calculate_planck_rosseland_opacity_block)
        call blk_traversal(calculate_fld_mfp_sigma_block)
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            if (resolve_rad_dt) blk%Erad_int=0d0
            blk=>blk%blk_next
        end do
        call calculate_fld_conductivity()
        dt=time_sys%dt_hydro
        if (resolve_rad_dt) then
            dtemp=50;iter=1
            do while (dt>0)
                call blk_traversal(get_cv_block)
                call petsc_fld_estimate_dt(dtemp)
                ddt=min(time_sys%dt_radiation,dt)
                call petsc_assemble_mat_vec(ddt,A,b)
                call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
                call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
                call VecAssemblyBegin(b,ierr);CHKERRQ(ierr)
                call VecAssemblyEnd(b,ierr);CHKERRQ(ierr)
                !call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
                !call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
                !call MatGetOwnershipRange(A,m,n,ierr)
                !print *,rank,m,n
                !stop
                call petsc_set_ksp(ksp,petsc_rtol)
                !call kspsetreusepreconditioner(ksp,PETSC_TRUE,ierr)

                call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)
                call petsc_fld_check_positivity(x,A)
                j=1;petsc_rtol_temp=petsc_rtol
                do while (processor%petsc_decrease_rtol)
                    if (j>5) then
                        print *,'petsc does not converge1'
                        stop
                    end if
                    petsc_rtol_temp=petsc_rtol_temp/10
                    call petsc_set_ksp(ksp,petsc_rtol_temp)
                    call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)
                    call petsc_fld_check_positivity(x)
                    j=j+1
                end do
                call petsc_accept_result(x)
                blk=>blk_processor_head
                do i=1,np_nblk(rank+1)
                    blk%Erad_int=blk%Erad_int+blk%Erad*ddt/time_sys%dt_hydro
                    if (i/=np_nblk(rank+1)) then
                        blk=>blk%blk_next
                    end if
                end do
                dt=dt-ddt
                iter=iter+1
                if (iter==20000) then
                    if (rank==0) print *, 'explicit subcycle, converge slowly'
                    stop
                end if
            end do
        else
#if         ieos==2
            call petsc_fld_predict_E_T(dt,A,b,ksp,x)
            call petsc_calculate_beta(dt)
            call petsc_fld_update_linear_system(A,b)
            call petsc_accept_result(x)
#else
            call blk_traversal(get_cv_block)
            call petsc_assemble_mat_vec(dt,A,b)
                
            call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
            call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
            call VecAssemblyBegin(b,ierr);CHKERRQ(ierr)
            call VecAssemblyEnd(b,ierr);CHKERRQ(ierr)
            call petsc_set_ksp(ksp,petsc_rtol)

            call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)
            call petsc_fld_check_positivity(x)
            j=1;petsc_rtol_temp=petsc_rtol
            do while (processor%petsc_decrease_rtol)
                if (j>5) then
                    print *,'petsc does not converge1'
                    stop
                end if
                petsc_rtol_temp=petsc_rtol_temp/10
                call petsc_set_ksp(ksp,petsc_rtol_temp)
                call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)
                call petsc_fld_check_positivity(x)
                j=j+1
            end do
            call petsc_accept_result(x)
#endif
        end if
        nullify(blk)
        call blk_traversal(convert_rho_temp_to_u_w_block)
        call communicate_blocks_fld_1d()
        call blk_traversal(petsc_fld_momentum)
        !call communicate_hydro()
    else if (nd==2) then
        dt=time_sys%dt_hydro
        call blk_traversal(get_cv_block)
        call blk_traversal(calculate_planck_rosseland_opacity_block)
        call blk_traversal(calculate_fld_mfp_sigma_block)
        call calculate_fld_conductivity()
        call petsc_assemble_mat_vec_2d(dt,A,b)

        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
        call VecAssemblyBegin(b,ierr);CHKERRQ(ierr)
        call VecAssemblyEnd(b,ierr);CHKERRQ(ierr)
        !call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
        !call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call petsc_set_ksp(ksp,petsc_rtol)
        call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)
        !call petsc_fld_check_positivity(x)
        call petsc_accept_result(x)
        call blk_traversal(convert_rho_temp_to_u_w_block)
        call communicate_blocks_2d_fld()
        !call petsc_fld_momentum(dt)
        !call communicate_hydro()
    end if
end subroutine petsc_fld_solve

subroutine petsc_fld_predict_E_T(dt,A,b,ksp,x)
    !solve the radiation transfer equation for once
    !save the solution of the radiation energy and gas temperature to E1 and T1
    type(blockdef), pointer :: blk
    type(axbsystem) userctx
    PetscScalar, pointer :: g(:)
    KSP  ksp
    Vec  b,x
    Mat  A
    real(8) :: dt,petsc_rtol_temp
    integer :: i,j,II,ierr
    call blk_traversal(get_cv_block)
    call petsc_assemble_mat_vec(dt,A,b)
        
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call VecAssemblyBegin(b,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(b,ierr);CHKERRQ(ierr)
    call petsc_set_ksp(ksp,petsc_rtol)

    call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)
    call petsc_fld_check_positivity(x)
    j=1;petsc_rtol_temp=petsc_rtol
    do while (processor%petsc_decrease_rtol)
        if (j>5) then
            print *,'petsc does not converge predict ET'
            stop
        end if
        petsc_rtol_temp=petsc_rtol_temp/10
        call petsc_set_ksp(ksp,petsc_rtol_temp)
        call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)
        call petsc_fld_check_positivity(x)
        j=j+1
    end do
    call VecGetArrayReadF90(x,g,ierr)
    do i=1,n_local_T
        II=2*i
        petsc_local_E1(i)=g(II-1)
        petsc_local_T1(i)=g(II)
    end do
    call VecRestoreArrayReadF90(x,g,ierr)
end subroutine petsc_fld_predict_E_T

subroutine petsc_calculate_beta(dt)
    type(blockdef), pointer :: blk
    real(8) :: dt
    integer :: i,j,ii
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            ii=(i-1)*blk_size_nx+j
            beta_rank(ii)=blk%sigma_planck(j,1,1)*c_light*dt
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end subroutine petsc_calculate_beta

subroutine petsc_fld_update_linear_system(A,b)
    !use E1 and T1 to update the linear system
    type(blockdef), pointer :: blk
    type(axbsystem) userctx
    Vec  b
    Mat  A
    real(8) :: v,rho,temp,egv,egvm,cv
    integer :: i,j,k,ii,jj,ierr
    do i=1,n_local_E
        ii=(i-1)*2
        jj=ii+1
        temp=petsc_local_T1(i)
        v=-4d0*beta_rank(i)*a_rad*temp**3
        call MatSetValue(A,ii,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=petsc_local_E1(i)-3d0*beta_rank(i)*a_rad*temp**4
        call VecSetValue(b,ii,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    end do
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            k=(i-1)*blk_size_nx+j
            ii=k*2-1
            egv=blk%egv(j,1,1)
            rho=blk%w(j,1,1,1)
            temp=petsc_local_T1(k)
#if         ieos==2
            cv=eos_h_cv(rho,temp)
            egvm=egvrhot(rho,temp)
#else
            cv=rho*kb/(gamma_gas-1)/maw/amu
            egvm=cv*temp
#endif
            v=cv+4d0*beta_rank(k)*a_rad*temp**3
            call MatSetValue(A,ii,ii,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
            v=cv*temp+3d0*beta_rank(k)*a_rad*temp**4+egv-egvm
            call VecSetValue(b,ii,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call VecAssemblyBegin(b,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(b,ierr);CHKERRQ(ierr)
end subroutine petsc_fld_update_linear_system

subroutine petsc_fld_check_positivity(x,A)
    type(blockdef), pointer :: blk
    Vec x
    Mat, optional :: A
    PetscScalar, pointer :: g(:)
    integer :: i,j,ierr,blk_id,II
    processor%processor_logical=.false.
    processor%global_logical=.false.
    call VecGetArrayReadF90(x,g,ierr)
    blk=>blk_processor_head
    outer:  do j=1,np_nblk(rank+1)
        blk_id=blk%blk_id-blk_processor_head%blk_id+1
        do i=1,blk_size_nx
            II=2*((blk_id-1)*blk_size_nx+i-1)+1
            if (g(II)<=0.or.g(II+1)<=0) then
                processor%processor_logical=.false.
                !print *,II,g(II),g(II+1),blk%temp(i,1,1)
                !if (g(II+1)<=0) then
                !    print *,'tgas<=0'
                !    print *,II,g(II),g(II+1),blk%temp(i,1,1)
                !    stop
                !else
                !    g(II)=1e-1!a_rad*g(II+1)**4
                !end if
                !processor%processor_logical=.true.
                exit outer
            else if (isnan(g(II)).or.g(II)>huge(g(II)).or.isnan(g(II+1)).or.g(II+1)>huge(g(II+1))) then
                print *,'Nan or Infinity in petsc',rank,j
                !print *,II,g(II),g(II+1),blk%temp(i,1,1)
                !processor%processor_logical=.false.
                !exit outer
                !call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)
                call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
                stop
            else
                processor%processor_logical=.true.
            end if
        end do
        blk=>blk%blk_next
    end do outer
    nullify(blk)
    call VecRestoreArrayReadF90(x,g,ierr)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_gather(processor%processor_logical,1,MPI_LOGICAL,processor%global_logical,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (rank==0) then
        processor%petsc_decrease_rtol=.not.(all(processor%global_logical))
    end if
    call mpi_bcast(processor%petsc_decrease_rtol,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
end subroutine petsc_fld_check_positivity

subroutine petsc_accept_result(x)
    type(blockdef), pointer :: blk
    Vec x
    PetscScalar, pointer :: g(:)
    integer :: i,j,k,ierr,blk_id,II
    real(8) :: maxtemp_pre,maxtemp_after,r1,r2,q,dx
    !maxtemp_pre=0;maxtemp_after=0
    call VecGetArrayReadF90(x,g,ierr)
    if (nd==1) then
        blk=>blk_processor_head
        do j=1,np_nblk(rank+1)
            blk_id=blk%blk_id-blk_processor_head%blk_id+1
            do i=1,blk_size_nx
                II=2*((blk_id-1)*blk_size_nx+i-1)+1
                blk%Erad(i,1,1)=g(II)
                blk%temp(i,1,1)=g(II+1)
            end do
            if (associated(blk,blk_head)) then
                blk%Erad(0,1,1)=blk%Erad(1,1,1)
            end if
            if (associated(blk,blk_tail)) then
                dx=blk%dxyz(1)
                if (igeometry==0) then
                    if (rad_bound_type(2)==2) then
                        blk%Erad(0,1,1)=blk%Erad(1,1,1)
                    else if (rad_bound_type(2)==1) then
                        q=1.5d0*blk%sigma_rosseland(blk_size_nx,1,1)*dx
                        blk%Erad(0,1,1)=blk%Erad(1,1,1)/(1d0+q)
                    end if
                else if (igeometry==2) then
                    if (rad_bound_type(2)==2) then
                        blk%Erad(blk_size_nx+1,1,1)=blk%Erad(blk_size_nx,1,1)
                    else if (rad_bound_type(2)==1) then
                        r1=blk%x_center(blk_size_nx)
                        r2=blk%x_center(blk_size_nx+1)
                        blk%Erad(blk_size_nx+1,1,1)=r1**2/r2**2*blk%Erad(blk_size_nx,1,1)
                    end if
                end if
            end if
            if (j/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
    else if (nd==2) then
        blk=>blk_processor_head
        do k=1,np_nblk(rank+1)
            blk_id=blk%blk_id-blk_processor_head%blk_id+1
            do j=1,blk_size_ny
                do i=1,blk_size_nx
                    II=2*((blk_id-1)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)+1
                    blk%Erad(i,j,1)=g(II)
                    blk%temp(i,j,1)=g(II+1)
                end do
            end do
            !print *,blk%Erad
            !print *,blk%temp
            if (k/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
    end if
    nullify(blk)
    call VecRestoreArrayReadF90(x,g,ierr)
end subroutine petsc_accept_result

subroutine petsc_set_ksp(ksp,rtol)
    KSP ksp
    KSP, allocatable, dimension(:) :: subksp
    PC pc,subpc
    real(8) :: rtol
    integer :: ierr,i,nlocal,first
    if (np==1) then
        call KSPSetType(ksp,KSPGMRES,ierr);CHKERRQ(ierr)
        call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
        call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)
        call PCSetType(pc,PCLU,ierr);CHKERRQ(ierr)
    else
        call KSPSetType(ksp,KSPGMRES,ierr);CHKERRQ(ierr)
        call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
        call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)
        call PCSetType(pc,PCASM,ierr)
        call PCASMSetType(pc,PC_ASM_INTERPOLATE,ierr)
        call PCASMSetOverlap(pc,1,ierr)
        call KSPsetup(ksp,ierr)
        call PCASMGetSubKSP(pc,nlocal,first,PETSC_NULL_KSP,ierr)
        allocate(subksp(nlocal))
        call PCASMGetSubKSP(pc,nlocal,first,subksp,ierr)
        do i=0,nlocal-1
            call KSPGetPC(subksp(i+1),subpc,ierr)
            call PCSetType(subpc,PCLU,ierr); CHKERRA(ierr)
            call KSPSetType(subksp(i+1),KSPGMRES,ierr); CHKERRA(ierr)
            call KSPSetTolerances(subksp(i+1),rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
        end do
        deallocate(subksp)
    end if
end subroutine petsc_set_ksp

subroutine petsc_assemble_mat_vec(dt,A,b)
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: dt,beta,dl,dr,cv,er,temp,alphal,alphar,gl,gr,dx,xl,xc,xr,r1,r2,   &
        r_mean,Sadv,eradl,eradr,frhol,frhor,rhol,rhor
    integer :: i,j,ii_E,ii_T,iblk
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        iblk=blk%blk_id
        dx=blk%dxyz(1)
        !first cell of this block
        if (associated(blk,blk_head)) then
            call petsc_Ab_1d_lb(dt,A,b)
        else
            ii_E=2*((iblk-1)*blk_size_nx+1-1)
            ii_T=ii_E+1
            beta=blk%sigma_planck(1,1,1)*c_light*dt
            dl=blk%kx(0,1,1)
            dr=blk%kx(1,1,1)
            cv=blk%cv(1,1,1)
            er=blk%Erad(1,1,1)
            temp=blk%temp(1,1,1)
            if (igeometry==0) then
                xl=blk%nb_coor%xl(1)
                xc=blk%x_center(1)
                alphal=dt/dx/(xc-xl)
                alphar=dt/dx/dx
                gl=1d0
                gr=1d0
            else if (igeometry==2) then
                xl=blk%nb_coor%xl(1)
                xc=blk%x_center(1)
                xr=blk%x_center(2)
                alphal=dt/dx/(xc-xl)
                alphar=dt/dx/(xr-xc)
                r1=blk%x_interface(0)
                r2=blk%x_interface(1)
                gl=r1**2/xc**2
                gr=r2**2/xc**2
            end if
            call petsc_Ab_1d_interior(beta,alphal,alphar,gl,gr,dl,dr,cv,er,temp,ii_E,ii_T,A,b)
        end if
        !interior cells of this block
        do j=2,blk_size_nx-1
            ii_E=2*((iblk-1)*blk_size_nx+j-1)
            ii_T=ii_E+1
            beta=blk%sigma_planck(j,1,1)*c_light*dt
            dl=blk%kx(j-1,1,1)
            dr=blk%kx(j,1,1)
            cv=blk%cv(j,1,1)
            er=blk%Erad(j,1,1)
            temp=blk%temp(j,1,1)
            if (igeometry==0) then
                alphal=dt/dx/dx
                alphar=dt/dx/dx
                gl=1d0
                gr=1d0
            else if (igeometry==2) then
                xl=blk%x_center(j-1)
                xc=blk%x_center(j)
                xr=blk%x_center(j+1)
                alphal=dt/dx/(xc-xl)
                alphar=dt/dx/(xr-xc)
                r1=blk%x_interface(j-1)
                r2=blk%x_interface(j)
                gl=r1**2/xc**2
                gr=r2**2/xc**2
            end if
            call petsc_Ab_1d_interior(beta,alphal,alphar,gl,gr,dl,dr,cv,er,temp,ii_E,ii_T,A,b)
        end do
        !last cell of this block
        if (associated(blk,blk_tail)) then
            call petsc_Ab_1d_rb(dt,A,b)
        else
            ii_E=2*((iblk-1)*blk_size_nx+blk_size_nx-1)
            ii_T=ii_E+1
            beta=blk%sigma_planck(blk_size_nx,1,1)*c_light*dt
            dl=blk%kx(blk_size_nx-1,1,1)
            dr=blk%kx(blk_size_nx,1,1)
            cv=blk%cv(blk_size_nx,1,1)
            er=blk%Erad(blk_size_nx,1,1)
            temp=blk%temp(blk_size_nx,1,1)
            if (igeometry==0) then
                xc=blk%x_center(blk_size_nx)
                xr=blk%nb_coor%xu(1)
                alphal=dt/dx/dx
                alphar=dt/dx/(xr-xc)
                gl=1d0
                gr=1d0
            else if (igeometry==2) then
                xl=blk%x_center(blk_size_nx-1)
                xc=blk%x_center(blk_size_nx)
                xr=blk%nb_coor%xu(1)
                alphal=dt/dx/(xc-xl)
                alphar=dt/dx/(xr-xc)
                r1=blk%x_interface(blk_size_nx-1)
                r2=blk%x_interface(blk_size_nx)
                gl=r1**2/xc**2
                gr=r2**2/xc**2
            end if
            call petsc_Ab_1d_interior(beta,alphal,alphar,gl,gr,dl,dr,cv,er,temp,ii_E,ii_T,A,b)
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    nullify(blk)
end subroutine petsc_assemble_mat_vec

subroutine petsc_Ab_1d_interior(beta,alphal,alphar,gl,gr,dl,dr,cv,er,temp,ii_E,ii_T,A,b)
    real(8) :: beta,alphal,alphar,gl,gr,dl,dr,cv,er,temp,v,sadv
    integer :: ii_E,ii_T,jj,ierr
    Mat A
    Vec b
    jj=ii_E-2
    v=-alphal*gl*dl*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=ii_E
    v=1+(alphal*gl*dl+alphar*gr*dr+beta)*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=ii_E+2
    v=-alphar*gr*dr*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4*beta*a_rad*temp**3+cv
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=er-(3*beta*a_rad*temp**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=cv*temp+3*beta*a_rad*temp**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_1d_interior

subroutine petsc_Ab_1d_lb(dt,A,b)
    type(blockdef), pointer :: blk
    real(8) :: dt,beta,alphal,alphar,gl,gr,dl,dr,cv,er,temp,v,dx,xc,xr,r1,Sadv
    integer :: ii_E,ii_T,jj,ierr
    Mat A
    Vec b
    !only zero gradient left boundary is implemented
    if (rank==0) then
        blk=>blk_head
        ii_E=0
        ii_T=1
        beta=blk%sigma_planck(1,1,1)*c_light*dt
        dr=blk%kx(1,1,1)
        cv=blk%cv(1,1,1)
        er=blk%Erad(1,1,1)
        temp=blk%temp(1,1,1)
        dx=blk%dxyz(1)
        if (igeometry==0) then
            alphar=dt/dx/dx
            gr=1d0
        else if (igeometry==2) then
            xc=blk%x_center(1)
            xr=blk%x_center(2)
            alphar=dt/dx/(xr-xc)
            r1=blk%x_interface(1)
            gr=r1**2/xc**2
        end if
        jj=ii_E
        v=1+(alphar*gr*dr+beta)*const_Erad
        call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        jj=ii_E+2
        v=-alphar*gr*dr*const_Erad
        call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=-4*beta*a_rad*temp**3*const_Erad
        call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=4*beta*a_rad*temp**3+cv
        call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=-beta
        call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=er-(3*beta*a_rad*temp**4)*const_Erad
        call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=cv*temp+3*beta*a_rad*temp**4
        call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        nullify(blk)
    end if
end subroutine petsc_Ab_1d_lb

subroutine petsc_Ab_1d_rb(dt,A,b)
    type(blockdef), pointer :: blk
    real(8) :: dt,beta,alphal,alphar,gl,gr,dl,dr,cv,er,temp,v,dx,q,xl,xc,xr,r1,r2,Sadv
    integer :: ii_E,ii_T,jj,ierr
    Mat A
    Vec b
    if (rank==np-1) then
        blk=>blk_tail
        ii_E=2*((nblk_total-1)*blk_size_nx+blk_size_nx-1)
        ii_T=ii_E+1
        beta=blk%sigma_planck(blk_size_nx,1,1)*c_light*dt
        dl=blk%kx(blk_size_nx-1,1,1)
        dr=blk%kx(blk_size_nx,1,1)
        cv=blk%cv(blk_size_nx,1,1)
        er=blk%Erad(blk_size_nx,1,1)
        temp=blk%temp(blk_size_nx,1,1)
        dx=blk%dxyz(1)
        if (igeometry==0) then
            alphal=dt/dx/dx
            alphar=dt/dx/dx
            gl=1d0
            gr=1d0
            if (rad_bound_type(2)==1) then
                q=1.5d0*blk%sigma_rosseland(blk_size_nx,1,1)*dx
                v=1+(beta+alphal*dl+q*alphar*dr/(1d0+q))*const_Erad
            else if (rad_bound_type(2)==2) then
                v=1+(beta+alphal*dl)*const_Erad
            end if
        else if (igeometry==2) then
            xl=blk%x_center(blk_size_nx-1)
            xc=blk%x_center(blk_size_nx)
            xr=blk%x_center(blk_size_nx+1)
            alphal=dt/dx/(xc-xl)
            alphar=dt/dx/(xr-xc)
            r1=blk%x_interface(blk_size_nx-1)
            r2=blk%x_interface(blk_size_nx)
            gl=r1**2/xc**2
            gr=r2**2/xc**2
            if (rad_bound_type(2)==1) then
                v=1+(alphal*gl*dl+alphar*gr*dr*(xr**2-xc**2)/xr**2+beta)
            else if (rad_bound_type(2)==2) then
                v=1+(alphal*gl*dl+beta)*const_Erad
            end if
        end if
        jj=ii_E
        call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        jj=ii_E-2
        v=-alphal*gl*dl*const_Erad
        call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=-4*beta*a_rad*temp**3*const_Erad
        call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=4*beta*a_rad*temp**3+cv
        call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=-beta
        call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=er-(3*beta*a_rad*temp**4)*const_Erad
        call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        v=cv*temp+3*beta*a_rad*temp**4
        call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        nullify(blk)
    end if
end subroutine petsc_Ab_1d_rb

subroutine petsc_fld_momentum(blk)
    type(blockdef), pointer :: blk
    real(8) :: dx,gradEx,Eext,Favg,Fdiff,x1,x2
    integer :: i,j,ijk(3)
    !the left boundary interface of this block
    if (blk%loc_type==1) then
        if (rad_bound_type(1)==2) then
            blk%Fradx(0,1,1)=0d0
        end if
    else
        x1=blk%nb_coor%xl(1)
        x2=blk%x_center(1)
        if (resolve_rad_dt) then
            gradEx=(blk%Erad_int(1,1,1)-blk%Erad_pre_int)/(x2-x1)
        else
            gradEx=(blk%Erad(1,1,1)-blk%Erad_pre)/(x2-x1)
        end if
        blk%Fradx(0,1,1)=-blk%kx(0,1,1)*gradEx
    end if
    !the right boundary interface of this block
    if (blk%loc_type==3) then
        x1=blk%x_center(blk_size_nx)
        x2=blk%x_center(blk_size_nx+1)
        if (resolve_rad_dt) then
            gradEx=(blk%Erad_int(blk_size_nx+1,1,1)-blk%Erad_int(blk_size_nx,1,1))/(x2-x1)
        else
            gradEx=(blk%Erad(blk_size_nx+1,1,1)-blk%Erad(blk_size_nx,1,1))/(x2-x1)
        end if
        blk%Fradx(blk_size_nx,1,1)=-blk%kx(blk_size_nx,1,1)*gradEx
    else
        x1=blk%x_center(blk_size_nx)
        x2=blk%nb_coor%xu(1)
        if (resolve_rad_dt) then
            gradEx=(blk%Erad_next_int-blk%Erad_int(blk_size_nx,1,1))/(x2-x1)
        else
            gradEx=(blk%Erad_next-blk%Erad(blk_size_nx,1,1))/(x2-x1)
        end if
        blk%Fradx(blk_size_nx,1,1)=-blk%kx(blk_size_nx,1,1)*gradEx
    end if
    !the interior interfaces of this block
    do i=1,blk_size_nx-1
        x1=blk%x_center(i)
        x2=blk%x_center(i+1)
        if (resolve_rad_dt) then
            gradEx=(blk%Erad_int(i+1,1,1)-blk%Erad_int(i,1,1))/(x2-x1)
        else
            gradEx=(blk%Erad(i+1,1,1)-blk%Erad(i,1,1))/(x2-x1)
        end if
        blk%Fradx(i,1,1)=-blk%kx(i,1,1)*gradEx
    end do
end subroutine petsc_fld_momentum

subroutine petsc_frad_consistency()
    type(blockdef), pointer :: blk,blk_next
    real(8) :: frad1,frad2,kx1,kx2,Erad1,Erad2,sigma1,sigma2
    if (nd==1) then
        blk=>blk_processor_head
        print *,'head',rank,blk%Fradx(0,1,1)
        blk=>blk_processor_tail
        print *,'tail',rank,blk%Fradx(blk_size_nx,1,1)
    else if (nd==2) then
    end if
end subroutine petsc_frad_consistency

subroutine petsc_fld_recycle(userctx,ierr)
    PetscErrorCode ierr
    type(axbsystem) userctx
    call VecDestroy(userctx%x,ierr);CHKERRQ(ierr)
    call VecDestroy(userctx%b,ierr);CHKERRQ(ierr)
    call MatDestroy(userctx%A,ierr);CHKERRQ(ierr)
    call KSPDestroy(userctx%ksp,ierr);CHKERRQ(ierr)
#if     ieos==2
    deallocate(petsc_local_T1,petsc_local_T2,petsc_local_E1,petsc_local_E2,petsc_diffT,petsc_diffE,beta_rank)
#endif
end subroutine petsc_fld_recycle

subroutine petsc_fld_finalize(userctx,ierr)
    PetscErrorCode ierr
    type(axbsystem) userctx
    if (refine_type=='static'.or.refine_type=='none') then
        call VecDestroy(userctx%x,ierr);CHKERRQ(ierr)
        call VecDestroy(userctx%b,ierr);CHKERRQ(ierr)
        call MatDestroy(userctx%A,ierr);CHKERRQ(ierr)
        call KSPDestroy(userctx%ksp,ierr);CHKERRQ(ierr)
    else if (refine_type=='adaptive'.or.refine_type=='mixed') then
    end if
end subroutine petsc_fld_finalize

end module petsc_fld
