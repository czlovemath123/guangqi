module radiation_common_functions
use phylib
use mathlib
use communication
use datastructure
use eos
implicit none

contains

subroutine generate_rosseland_planck_kap_table()
    call readin_gas_opacity_tables()
    call readin_dust_opacity()
end subroutine generate_rosseland_planck_kap_table

subroutine readin_gas_opacity_tables()
    real(8), dimension(:,:,:), allocatable :: kappa_rosseland,kappa_planck
    real(8), dimension(:,:), allocatable :: kappa_R_trim,kappa_P_trim
    real(8), dimension(:), allocatable :: tgas_table,pres_table
    real(8) :: metal,pres,rho,kR,kP,metal_list(3),tgas_list(94),pres_list(126)
    character(len=64) :: fmt1
    integer :: i,j,k,nt,npres,nmetal,tgas,trad,ipres_min,ipres_max,it_min,it_max,nt_table,npres_table
    nt=94
    npres=126
    nmetal=3
    fmt1='(1F5.1,1I8,2E10.4,2E13.6)'
    opacity_gas_pres_min=calculate_p_from_rho_t(opacity_gas_rho_min,opacity_gas_t_min)
    opacity_gas_pres_max=calculate_p_from_rho_t(opacity_gas_rho_max,opacity_gas_t_max)
    allocate(kappa_rosseland(nmetal,nt,npres),kappa_planck(nmetal,nt,npres))
    open(unit=11,file='tables/opacity/table1T.dat',status='old',action='read')
    do i=1,nmetal
        do j=1,nt
            do k=1,npres
                read(11,fmt=fmt1) metal,tgas,pres,rho,kR,kP
                kappa_rosseland(i,j,k)=kR
                kappa_planck(i,j,k)=kP
                if (i==metallicity.and.j==1) then
                    pres_list(k)=pres
                end if
            end do
            if (i==metallicity) then
                tgas_list(j)=tgas
            end if
        end do
        metal_list(i)=metal
    end do
    close(11)
    do i=1,npres-1
        if (pres_list(i+1)>opacity_gas_pres_min) then
            ipres_min=i
            exit
        end if
    end do
    do i=npres,2,-1
        if (pres_list(i-1)<opacity_gas_pres_max) then
            ipres_max=i
            exit
        end if
    end do
    npres_table=ipres_max-ipres_min+1
    do i=1,nt-1
        if (tgas_list(i+1)>opacity_gas_t_min) then
            it_min=i
            exit
        end if
    end do
    do i=nt,2,-1
        if (tgas_list(i-1)<opacity_gas_t_max) then
            it_max=i
            exit
        end if
    end do
    nt_table=it_max-it_min+1
    call create_table_2d(nt_table,npres_table,rosseland_gas_opacity_table)
    call create_table_2d(nt_table,npres_table,planck_gas_opacity_table)
    allocate(pres_table(npres_table),tgas_table(nt_table))
    allocate(kappa_R_trim(nt_table,npres_table),kappa_P_trim(nt_table,npres_table))
    do i=ipres_min,ipres_max
        pres_table(i-ipres_min+1)=pres_list(i)
    end do
    do i=it_min,it_max
        tgas_table(i-it_min+1)=tgas_list(i)
    end do
    do i=ipres_min,ipres_max
        do j=it_min,it_max
            kappa_R_trim(j-it_min+1,i-ipres_min+1)=kappa_rosseland(2,j,i)
            kappa_P_trim(j-it_min+1,i-ipres_min+1)=kappa_planck(2,j,i)
        end do
    end do
    kappa_P_trim=log10(kappa_P_trim)
    kappa_R_trim=log10(kappa_R_trim)
    pres_table=log10(pres_table)
    tgas_table=log10(tgas_table)
    rosseland_gas_opacity_table%t2d=kappa_R_trim
    rosseland_gas_opacity_table%xlist=tgas_table
    rosseland_gas_opacity_table%ylist=pres_table
    rosseland_gas_opacity_table%nd1=nt_table
    rosseland_gas_opacity_table%nd2=npres_table
    planck_gas_opacity_table%t2d=kappa_P_trim
    planck_gas_opacity_table%xlist=tgas_table
    planck_gas_opacity_table%ylist=pres_table
    planck_gas_opacity_table%nd1=nt_table
    planck_gas_opacity_table%nd1=npres_table
    deallocate(pres_table,tgas_table,kappa_R_trim,kappa_P_trim)
end subroutine readin_gas_opacity_tables

subroutine readin_gas_opacity_two_temp_tables()
end subroutine readin_gas_opacity_two_temp_tables

subroutine readin_dust_opacity()
    real(8), dimension(:,:), allocatable :: kappa_rosseland,kappa_planck
    real(8) :: rho,t,kp,kr,t_table(101),rho_table(61)
    integer :: nt,nrho,i,j,k
    character(len=64) :: fmt1
    nt=101
    nrho=61
    fmt1='(3ES18.5E2)'
    allocate(kappa_rosseland(nt,nrho),kappa_planck(nt,nrho))
    open(unit=11,file='tables/opacity/kP.dat',status='old',action='read')
    open(unit=12,file='tables/opacity/kR.dat',status='old',action='read')
    do i=1,nrho
        do j=1,nt
            read(11,fmt=fmt1) rho,t,kp
            read(12,fmt=fmt1) rho,t,kr
            if (i==1) then
                t_table(j)=t
            end if
            kappa_rosseland(j,i)=kr
            kappa_planck(j,i)=kp
        end do
        rho_table(i)=rho
    end do
    close(11)
    close(12)
    call create_table_2d(nt,nrho,rosseland_dust_opacity_table)
    call create_table_2d(nt,nrho,planck_dust_opacity_table)
    t_table=log10(t_table)
    rho_table=log10(rho_table)
    kappa_rosseland=log10(kappa_rosseland)
    rosseland_dust_opacity_table%t2d=kappa_rosseland
    rosseland_dust_opacity_table%xlist=t_table
    rosseland_dust_opacity_table%ylist=rho_table
    rosseland_dust_opacity_table%nd1=nt
    rosseland_dust_opacity_table%nd2=nrho
    kappa_planck=log10(kappa_planck)
    planck_dust_opacity_table%t2d=kappa_planck
    planck_dust_opacity_table%xlist=t_table
    planck_dust_opacity_table%ylist=rho_table
    planck_dust_opacity_table%nd1=nt
    planck_dust_opacity_table%nd1=nrho
end subroutine readin_dust_opacity

subroutine generate_rosseland_planck_kap_table_mg()
    call readin_gas_opacity_tables_mg()
    call readin_dust_opacity_mg()
end subroutine generate_rosseland_planck_kap_table_mg

subroutine readin_gas_opacity_tables_mg()
    real(8), dimension(:,:,:), allocatable :: kappa_rosseland,kappa_planck
    real(8), dimension(:,:), allocatable :: kappa_R_trim,kappa_P_trim
    real(8), dimension(:), allocatable :: tgas_table,rho_table
    real(8) :: metal,pres,rho,kR,kP,metal_list(3),tgas_list(94),rho_list(126)
    character(len=64) :: fmt1
    integer :: i,j,k,nt,nrho,nmetal,tgas,trad,irho_min,irho_max,it_min,it_max,nt_table,nrho_table
    nt=94
    nrho=126
    nmetal=3
    fmt1='(1F5.1,1I8,2E10.4,2E13.6)'
    allocate(kappa_rosseland(nmetal,nt,nrho),kappa_planck(nmetal,nt,nrho))
    open(unit=11,file='tables/opacity/table1T.dat',status='old',action='read')
    do i=1,nmetal
        do j=1,nt
            do k=1,nrho
                read(11,fmt=fmt1) metal,tgas,pres,rho,kR,kP
                kappa_rosseland(i,j,k)=kR
                kappa_planck(i,j,k)=kP
                if (i==2.and.j==1) then
                    rho_list(k)=rho
                end if
            end do
            if (i==2) then
                tgas_list(j)=tgas
            end if
        end do
        metal_list(i)=metal
    end do
    close(11)
    do i=1,nrho-1
        if (rho_list(i+1)>opacity_gas_rho_min) then
            irho_min=i
            exit
        end if
    end do
    do i=nrho,2,-1
        if (rho_list(i-1)<opacity_gas_rho_max) then
            irho_max=i
            exit
        end if
    end do
    nrho_table=irho_max-irho_min+1
    do i=1,nt-1
        if (tgas_list(i+1)>opacity_gas_t_min) then
            it_min=i
            exit
        end if
    end do
    do i=nt,2,-1
        if (tgas_list(i-1)<opacity_gas_t_max) then
            it_max=i
            exit
        end if
    end do
    nt_table=it_max-it_min+1
    call create_table_2d(nt_table,nrho_table,rosseland_gas_opacity_table)
    call create_table_2d(nt_table,nrho_table,planck_gas_opacity_table)
    allocate(rho_table(nrho_table),tgas_table(nt_table))
    allocate(kappa_R_trim(nt_table,nrho_table),kappa_P_trim(nt_table,nrho_table))
    do i=irho_min,irho_max
        rho_table(i-irho_min+1)=rho_list(i)
    end do
    do i=it_min,it_max
        tgas_table(i-it_min+1)=tgas_list(i)
    end do
    do i=irho_min,irho_max
        do j=it_min,it_max
            kappa_R_trim(j-it_min+1,i-irho_min+1)=kappa_rosseland(2,j,i)
            kappa_P_trim(j-it_min+1,i-irho_min+1)=kappa_planck(2,j,i)
        end do
    end do
    rho_table=log10(rho_table)
    tgas_table=log10(tgas_table)
    rosseland_gas_opacity_table%t2d=kappa_R_trim
    rosseland_gas_opacity_table%xlist=tgas_table
    rosseland_gas_opacity_table%ylist=rho_table
    rosseland_gas_opacity_table%nd1=nt_table
    rosseland_gas_opacity_table%nd2=nrho_table
    planck_gas_opacity_table%t2d=kappa_P_trim
    planck_gas_opacity_table%xlist=tgas_table
    planck_gas_opacity_table%ylist=rho_table
    planck_gas_opacity_table%nd1=nt_table
    planck_gas_opacity_table%nd1=nrho_table
    deallocate(rho_table,tgas_table,kappa_R_trim,kappa_P_trim)
end subroutine readin_gas_opacity_tables_mg

subroutine readin_dust_opacity_mg()
    real(8), dimension(:,:), allocatable :: kappa_rosseland,kappa_planck
    real(8) :: rho,t,kp,kr,t_table(101),rho_table(61)
    integer :: nt,nrho,i,j,k
    character(len=64) :: fmt1
    nt=101
    nrho=61
    fmt1='(3ES18.5E2)'
    allocate(kappa_rosseland(nt,nrho),kappa_planck(nt,nrho))
    open(unit=11,file='tables/opacity/kP.dat',status='old',action='read')
    open(unit=12,file='tables/opacity/kR.dat',status='old',action='read')
    do i=1,nrho
        do j=1,nt
            read(11,fmt=fmt1) rho,t,kp
            read(12,fmt=fmt1) rho,t,kr
            if (i==1) then
                t_table(j)=t
            end if
            kappa_rosseland(j,i)=kr
            kappa_planck(j,i)=kp
        end do
        rho_table(i)=rho
    end do
    close(11)
    close(12)
    call create_table_2d(nt,nrho,rosseland_dust_opacity_table)
    call create_table_2d(nt,nrho,planck_dust_opacity_table)
    t_table=log10(t_table)
    rho_table=log10(rho_table)
    rosseland_dust_opacity_table%t2d=kappa_rosseland
    rosseland_dust_opacity_table%xlist=t_table
    rosseland_dust_opacity_table%ylist=rho_table
    rosseland_dust_opacity_table%nd1=nt
    rosseland_dust_opacity_table%nd2=nrho
    planck_dust_opacity_table%t2d=kappa_planck
    planck_dust_opacity_table%xlist=t_table
    planck_dust_opacity_table%ylist=rho_table
    planck_dust_opacity_table%nd1=nt
    planck_dust_opacity_table%nd1=nrho
end subroutine readin_dust_opacity_mg

function fld_rosseland_opacity_interp(rho,tgas)
    real(8) :: rho,pres,tgas,fld_rosseland_opacity_interp,logrho,logpres,logtgas,kappa1,kappa2,kappa
#if     iopacity==1
    if (tgas<=200) then
        fld_rosseland_opacity_interp=1.87151d0
    else if (tgas>=1d6) then
        fld_rosseland_opacity_interp=1.085812E+01 
    else
        logtgas=log10(tgas)
        pres=calculate_p_from_rho_t(rho,tgas)
        logpres=log10(pres)
        logrho=log10(rho)
        if (tgas<=1500d0) then
            !use dust table
            logrho=min(max(logrho,-18d0),-8d0)
            call interpolation2d_linear(logtgas,logrho,kappa1,rosseland_dust_opacity_table%xlist,rosseland_dust_opacity_table%ylist,  &
                rosseland_dust_opacity_table%t2d)
        else
            kappa1=-10d0
        end if
        if (tgas>=1100d0) then
            !use gas table
            logpres=max(logpres,-9d0)
            call interpolation2d_linear(logtgas,logpres,kappa2,rosseland_gas_opacity_table%xlist,rosseland_gas_opacity_table%ylist,  &
                rosseland_gas_opacity_table%t2d)
        else
            kappa2=-10d0
        end if
        kappa=max(kappa1,kappa2)
        fld_rosseland_opacity_interp=10d0**kappa
    end if
#elif   iopacity==2
    fld_rosseland_opacity_interp=const_opacity
#endif
end function fld_rosseland_opacity_interp

subroutine fld_rosseland_opacity_interp_mg(rho,tgas,kappa_rosseland_mg)
    real(8) :: rho,tgas,logrho,logtgas,kappa
    real(8), dimension(:), allocatable :: kappa_rosseland_mg
!#if     iopacity==1
!    if (tgas<=700) then
!        kappa_rosseland_mg=1.87151d0
!    else if (tgas>=1d6) then
!        kappa_rosseland_mg=1.085812E+01 
!    else
!        logtgas=log10(tgas)
!        logrho=log10(rho)
!        if (tgas<=1500d0) then
!            !use dust table
!            logrho=min(max(logrho,-18d0),-8d0)
!            call interpolation2d_linear(logtgas,logrho,kappa,rosseland_dust_opacity_table%xlist,rosseland_dust_opacity_table%ylist,  &
!                rosseland_dust_opacity_table%t2d)
!        else
!            !use gas table
!            logrho=max(logrho,-20d0)
!            call interpolation2d_linear(logtgas,logrho,kappa,rosseland_gas_opacity_table%xlist,rosseland_gas_opacity_table%ylist,  &
!                rosseland_gas_opacity_table%t2d)
!        end if
!        kappa_rosseland_mg=kappa
!    end if
!#elif   iopacity==2
!    kappa_rosseland_mg=const_opacity
!#endif
end subroutine fld_rosseland_opacity_interp_mg

function fld_planck_opacity_interp_trad(rho,Erad,tgas)
    real(8) :: rho,pres,Erad,trad,tgas,fld_planck_opacity_interp_trad,logrho,logpres,logtrad,kappa,kappa1,kappa2
#if     iopacity==1
    trad=(Erad/a_rad)**0.25d0
    if (trad<=200) then
        fld_planck_opacity_interp_trad=2.77749d0
    else if (trad>=1d6) then
        fld_planck_opacity_interp_trad=7.693802E+00
    else
        logtrad=log10(trad)
        pres=calculate_p_from_rho_t(rho,tgas)
        logpres=log10(pres)
        logrho=log10(rho)
        if (trad<1500d0) then
            !use dust table
            logrho=min(max(logrho,-18d0),-8d0)
            call interpolation2d_linear(logtrad,logrho,kappa1,planck_dust_opacity_table%xlist,planck_dust_opacity_table%ylist,  &
                planck_dust_opacity_table%t2d)
        else
            kappa1=-10d0
        end if
        if (trad>1100d0) then
            !use gas table
            logpres=max(logpres,-9d0)
            call interpolation2d_linear(logtrad,logpres,kappa2,planck_gas_opacity_table%xlist,planck_gas_opacity_table%ylist,  &
                planck_gas_opacity_table%t2d)
        else
            kappa2=-10d0
        end if
        kappa=max(kappa1,kappa2)
        fld_planck_opacity_interp_trad=10d0**kappa
    end if
    !fld_planck_opacity_interp_trad=min(kappa,1d2)
#elif   iopacity==2
    fld_planck_opacity_interp_trad=const_opacity
#endif
end function fld_planck_opacity_interp_trad

function fld_planck_opacity_interp_tgas(rho,tgas)
    real(8) :: rho,pres,tgas,fld_planck_opacity_interp_tgas,logrho,logpres,logtgas,kappa
#if     iopacity==1
    if (tgas<=700) then
        fld_planck_opacity_interp_tgas=2.77749d0
    else if (tgas>=1d6) then
        fld_planck_opacity_interp_tgas=7.693802E+00
    else
        logtgas=log10(tgas)
        pres=calculate_p_from_rho_t(rho,tgas)
        logpres=log10(pres)
        logrho=log10(rho)
        if (tgas<1500d0) then
            !use dust table
            logrho=min(max(logrho,-18d0),-8d0)
            call interpolation2d_linear(logtgas,logrho,kappa,planck_dust_opacity_table%xlist,planck_dust_opacity_table%ylist,  &
                planck_dust_opacity_table%t2d)
        else
            !use gas table
            logpres=max(logpres,-9d0)
            call interpolation2d_linear(logtgas,logpres,kappa,planck_gas_opacity_table%xlist,planck_gas_opacity_table%ylist,  &
                planck_gas_opacity_table%t2d)
        end if
        fld_planck_opacity_interp_tgas=10d0**kappa
    end if
#elif   iopacity==2
    fld_planck_opacity_interp_tgas=const_opacity
#endif
end function fld_planck_opacity_interp_tgas

subroutine fld_planck_opacity_interp_mg(rho,tgas,kappa_planck_mg)
    real(8) :: rho,tgas,logrho,logtgas,kappa
    real(8), dimension(:), allocatable :: kappa_planck_mg
!#if     iopacity==1
!    if (tgas<=700) then
!        kappa_planck_mg=2.77749d0
!    else if (tgas>=1d6) then
!        kappa_planck_mg=7.693802E+00
!    else
!        logtgas=log10(tgas)
!        logrho=log10(rho)
!        if (tgas<1500d0) then
!            !use dust table
!            logrho=min(max(logrho,-18d0),-8d0)
!            call interpolation2d_linear(logtgas,logrho,kappa,planck_dust_opacity_table%xlist,planck_dust_opacity_table%ylist,  &
!                planck_dust_opacity_table%t2d)
!        else
!            !use gas table
!            logrho=max(logrho,-20d0)
!            call interpolation2d_linear(logtgas,logrho,kappa,planck_gas_opacity_table%xlist,planck_gas_opacity_table%ylist,  &
!                planck_gas_opacity_table%t2d)
!        end if
!        kappa_planck_mg=kappa
!    end if
!#elif   iopacity==2
!    kappa_planck_mg=const_opacity
!#endif
end subroutine fld_planck_opacity_interp_mg

subroutine generate_fld_mg()
    !group_frequency gives the lower and upper bound of each rad group
    !each group is further divided into nradcells
    integer :: i,j
    real(8), dimension(:), allocatable :: group_frequency
    real(8) :: dlnnu,multiply,tmin,tmax,numin,numax,dlnnu_group
    tmin=1d2
    tmax=1d5
    numin=planck_function_peak_frequency(tmin)
    numax=planck_function_peak_frequency(tmax)
    numin=6d13
    numax=3d15
    dlnnu_group=log(numax/numin)/nmg
    allocate(group_frequency(nmg+1))
    group_frequency(1)=numin
    do i=2,nmg+1
        group_frequency(i)=group_frequency(i-1)*exp(dlnnu_group)
    end do
    allocate(rad_mg(nradcell+1,nmg))
    do j=1,nmg
        dlnnu=log(group_frequency(j+1)/group_frequency(j))/nradcell
        multiply=exp(dlnnu)
        rad_mg(1,j)=group_frequency(j)
        do i=2,nradcell+1
            rad_mg(i,j)=rad_mg(i-1,j)*multiply
        end do
    end do
    deallocate(group_frequency)
end subroutine generate_fld_mg

function fld_flux_limiter(r)
    !lambda
    real(8) :: fld_flux_limiter,r
    if (r>=0.and.r<=1.5d0) then
        fld_flux_limiter=2d0/(3d0+sqrt(9d0+12d0*r*r))
    else
        fld_flux_limiter=1d0/(1d0+r+sqrt(1d0+2d0*r))
    end if
end function fld_flux_limiter

!function fld_flux_limiter(r)
!    real(8) :: fld_flux_limiter,r
!    fld_flux_limiter=(2d0+r)/(6d0+3d0*r+r**2)
!end function fld_flux_limiter

subroutine calculate_planck_rosseland_opacity_block(blk)
    type(blockdef), pointer :: blk
    real(8) :: Erad_cell,tgas,rho,kappa_planck_trad,kappa_planck_tgas
    integer :: i,j
    if (nd==1) then
        do i=blk_xlb,blk_xub
            Erad_cell=blk%Erad(i,1,1)
            tgas=blk%temp(i,1,1)
            rho=blk%w(i,1,1,1)
            kappa_planck_trad=fld_planck_opacity_interp_trad(rho,Erad_cell,tgas)
            !kappa_planck_tgas=fld_planck_opacity_interp_tgas(rho,tgas)
            blk%kappa_planck(i,1,1)=kappa_planck_trad!sqrt(kappa_planck_trad*kappa_planck_tgas)
            blk%kappa_rosseland(i,1,1)=fld_rosseland_opacity_interp(rho,tgas)
        end do
        if (associated(blk,blk_head)) then
            !left boundary block
            if (rad_bound_type(1)==9) then
                Erad_cell=blk%Erad(0,1,1)
                tgas=blk%temp(0,1,1)
                rho=blk%w(0,1,1,1)
                kappa_planck_trad=fld_planck_opacity_interp_trad(rho,Erad_cell,tgas)
                !kappa_planck_tgas=fld_planck_opacity_interp_tgas(rho,tgas)
                blk%kappa_planck(0,1,1)=kappa_planck_trad!sqrt(kappa_planck_trad*kappa_planck_tgas)
                blk%kappa_rosseland(0,1,1)=fld_rosseland_opacity_interp(rho,tgas)
            else
                blk%kappa_planck(0,1,1)=blk%kappa_planck(1,1,1)
                blk%kappa_rosseland(0,1,1)=blk%kappa_rosseland(1,1,1)
            end if
        end if
        if (associated(blk,blk_tail)) then
            !right boundary block
            if (rad_bound_type(2)==9) then
                Erad_cell=blk%Erad(blk_size_nx+1,1,1)
                tgas=blk%temp(blk_size_nx+1,1,1)
                rho=blk%w(blk_size_nx+1,1,1,1)
                kappa_planck_trad=fld_planck_opacity_interp_trad(rho,Erad_cell,tgas)
                !kappa_planck_tgas=fld_planck_opacity_interp_tgas(rho,tgas)
                blk%kappa_planck(blk_size_nx+1,1,1)=kappa_planck_trad!sqrt(kappa_planck_trad*kappa_planck_tgas)
                blk%kappa_rosseland(blk_size_nx+1,1,1)=fld_rosseland_opacity_interp(rho,tgas)
            else
                blk%kappa_planck(blk_size_nx+1,1,1)=blk%kappa_planck(blk_size_nx,1,1)
                blk%kappa_rosseland(blk_size_nx+1,1,1)=blk%kappa_rosseland(blk_size_nx,1,1)
            end if
        end if
#if     iopacity==2
        blk%kappa_planck=const_opacity
        blk%kappa_rosseland=const_opacity
#endif
    else if (nd==2) then
        do j=blk_ylb,blk_yub
            do i=blk_xlb,blk_xub
                Erad_cell=blk%Erad(i,j,1)
                tgas=blk%temp(i,j,1)
                rho=blk%w(i,j,1,1)
                !blk%kappa_planck(i,j,1)=fld_planck_opacity_interp_tgas(rho,tgas)
                blk%kappa_planck(i,j,1)=fld_planck_opacity_interp_trad(rho,Erad_cell,tgas)
                blk%kappa_rosseland(i,j,1)=fld_rosseland_opacity_interp(rho,tgas)
            end do
        end do
        select case (blk%loc_type)
        case (1)
            do i=1,blk_size_nx
                blk%kappa_planck(i,blk_size_ny+1,1)=blk%kappa_planck(i,blk_size_ny,1)
                blk%kappa_rosseland(i,blk_size_ny+1,1)=blk%kappa_rosseland(i,blk_size_ny,1)
            end do
            do i=1,blk_size_ny
                blk%kappa_planck(0,i,1)=blk%kappa_planck(1,i,1)
                blk%kappa_rosseland(0,i,1)=blk%kappa_rosseland(1,i,1)
            end do
        case (2)
            do i=1,blk_size_nx
                blk%kappa_planck(i,blk_size_ny+1,1)=blk%kappa_planck(i,blk_size_ny,1)
                blk%kappa_rosseland(i,blk_size_ny+1,1)=blk%kappa_rosseland(i,blk_size_ny,1)
            end do
        case (3)
            do i=1,blk_size_nx
                blk%kappa_planck(i,blk_size_ny+1,1)=blk%kappa_planck(i,blk_size_ny,1)
                blk%kappa_rosseland(i,blk_size_ny+1,1)=blk%kappa_rosseland(i,blk_size_ny,1)
            end do
            do i=1,blk_size_ny
                blk%kappa_planck(blk_size_nx+1,i,1)=blk%kappa_planck(blk_size_nx,i,1)
                blk%kappa_rosseland(blk_size_nx+1,i,1)=blk%kappa_rosseland(blk_size_nx,i,1)
            end do
        case (4)
            do i=1,blk_size_ny
                blk%kappa_planck(0,i,1)=blk%kappa_planck(1,i,1)
                blk%kappa_rosseland(0,i,1)=blk%kappa_rosseland(1,i,1)
            end do
        case (5)
            !do nothing
        case (6)
            do i=1,blk_size_ny
                blk%kappa_planck(blk_size_nx+1,i,1)=blk%kappa_planck(blk_size_nx,i,1)
                blk%kappa_rosseland(blk_size_nx+1,i,1)=blk%kappa_rosseland(blk_size_nx,i,1)
            end do
        case (7)
            do i=1,blk_size_nx
                blk%kappa_planck(i,0,1)=blk%kappa_planck(i,1,1)
                blk%kappa_rosseland(i,0,1)=blk%kappa_rosseland(i,1,1)
            end do
            do i=1,blk_size_ny
                blk%kappa_planck(0,i,1)=blk%kappa_planck(1,i,1)
                blk%kappa_rosseland(0,i,1)=blk%kappa_rosseland(1,i,1)
            end do
        case (8)
            do i=1,blk_size_nx
                blk%kappa_planck(i,0,1)=blk%kappa_planck(i,1,1)
                blk%kappa_rosseland(i,0,1)=blk%kappa_rosseland(i,1,1)
            end do
        case (9)
            do i=1,blk_size_nx
                blk%kappa_planck(i,0,1)=blk%kappa_planck(i,1,1)
                blk%kappa_rosseland(i,0,1)=blk%kappa_rosseland(i,1,1)
            end do
            do i=1,blk_size_ny
                blk%kappa_planck(blk_size_nx+1,i,1)=blk%kappa_planck(blk_size_nx,i,1)
                blk%kappa_rosseland(blk_size_nx+1,i,1)=blk%kappa_rosseland(blk_size_nx,i,1)
            end do
        end select
    end if
end subroutine calculate_planck_rosseland_opacity_block

subroutine calculate_planck_rosseland_opacity_mg_block(blk)
    type(blockdef), pointer :: blk
    real(8) :: tgas,rho
    real(8), dimension(:), allocatable :: kappa_planck_mg,kappa_rosseland_mg
    integer :: i,j
    if (nd==1) then
        do i=blk_xlb,blk_xub
            tgas=blk%temp(i,1,1)
            rho=blk%w(i,1,1,1)
            allocate(kappa_planck_mg(nmg),kappa_rosseland_mg(nmg))
            call fld_planck_opacity_interp_mg(rho,tgas,kappa_planck_mg)
            call fld_rosseland_opacity_interp_mg(rho,tgas,kappa_rosseland_mg)
            blk%kappa_planck_mg(:,i,1,1)=kappa_planck_mg
            blk%kappa_rosseland_mg(:,i,1,1)=kappa_rosseland_mg
            deallocate(kappa_planck_mg,kappa_rosseland_mg)
        end do
        if (associated(blk,blk_head)) then
            if (rad_bound_type(1)==9) then
                tgas=blk%temp(0,1,1)
                rho=blk%w(0,1,1,1)
                allocate(kappa_planck_mg(nmg),kappa_rosseland_mg(nmg))
                call fld_planck_opacity_interp_mg(rho,tgas,kappa_planck_mg)
                call fld_rosseland_opacity_interp_mg(rho,tgas,kappa_rosseland_mg)
                blk%kappa_planck_mg(:,0,1,1)=kappa_planck_mg
                blk%kappa_rosseland_mg(:,0,1,1)=kappa_rosseland_mg
                deallocate(kappa_planck_mg,kappa_rosseland_mg)
            else
                blk%kappa_planck_mg(:,0,1,1)=blk%kappa_planck_mg(:,1,1,1)
                blk%kappa_rosseland_mg(:,0,1,1)=blk%kappa_rosseland_mg(:,1,1,1)
            end if
        end if
        if (associated(blk,blk_tail)) then
            if (rad_bound_type(2)==9) then
                tgas=blk%temp(blk_size_nx+1,1,1)
                rho=blk%w(blk_size_nx+1,1,1,1)
                allocate(kappa_planck_mg(nmg),kappa_rosseland_mg(nmg))
                call fld_planck_opacity_interp_mg(rho,tgas,kappa_planck_mg)
                call fld_rosseland_opacity_interp_mg(rho,tgas,kappa_rosseland_mg)
                blk%kappa_planck_mg(:,blk_size_nx+1,1,1)=kappa_planck_mg
                blk%kappa_rosseland_mg(:,blk_size_nx+1,1,1)=kappa_rosseland_mg
                deallocate(kappa_planck_mg,kappa_rosseland_mg)
            else
                blk%kappa_planck_mg(:,blk_size_nx+1,1,1)=blk%kappa_planck_mg(:,blk_size_nx,1,1)
                blk%kappa_rosseland_mg(:,blk_size_nx+1,1,1)=blk%kappa_rosseland_mg(:,blk_size_nx,1,1)
            end if
        end if
#if     iopacity==2
        blk%kappa_planck_mg=const_opacity
        blk%kappa_rosseland_mg=const_opacity
#endif
    else if (nd==2) then
    end if
end subroutine calculate_planck_rosseland_opacity_mg_block

subroutine calculate_sc_opacity_block(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    do j=blk_ylb,blk_yub
        do i=blk_xlb,blk_xub
            blk%kappa_sc(i,j,1)=1d0
        end do
    end do
end subroutine calculate_sc_opacity_block

subroutine calculate_fld_mfp_sigma_block(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    if (nd==1) then
        do i=blk_xlb,blk_xub
            blk%sigma_planck(i,1,1)=blk%kappa_planck(i,1,1)*blk%w(i,1,1,1)
            blk%sigma_rosseland(i,1,1)=blk%kappa_rosseland(i,1,1)*blk%w(i,1,1,1)
        end do
        !blk%sigma_planck=3.1d-10
        !blk%sigma_rosseland=3.1d-10
        !blk%sigma_planck=0d0
        !blk%sigma_rosseland=1d7
    else if (nd==2) then
        !do j=blk_ylb,blk_yub
        !    do i=blk_xlb,blk_xub
        !        blk%sigma_planck(i,j,1)=blk%kappa_planck(i,j,1)*blk%w(i,j,1,1)
        !        blk%sigma_rosseland(i,j,1)=blk%kappa_rosseland(i,j,1)*blk%w(i,j,1,1)
        !    end do
        !end do
        !blk%sigma_planck=3.1d10
        !blk%sigma_rosseland=3.1d10
        blk%sigma_planck=1d2
        blk%sigma_rosseland=1d2
        !blk%sigma_planck=0d0
        !blk%sigma_rosseland=1d7
    end if
end subroutine calculate_fld_mfp_sigma_block

subroutine calculate_fld_mfp_sigma_mg_block(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    if (nd==1) then
        do i=blk_xlb,blk_xub
            do j=1,nmg
                blk%sigma_planck_mg(j,i,1,1)=blk%kappa_planck_mg(j,i,1,1)*blk%w(i,1,1,1)
                blk%sigma_rosseland_mg(j,i,1,1)=blk%kappa_rosseland_mg(j,i,1,1)*blk%w(i,1,1,1)
            end do
        end do
        blk%sigma_planck_mg=3.1d-10
        !blk%sigma_planck_mg(19:20,:,1,1)=0d0
        blk%sigma_rosseland_mg=3.1d-10
        !blk%sigma_planck_mg=0d0
        !blk%sigma_rosseland_mg=1d7
    else if (nd==2) then
    end if
end subroutine calculate_fld_mfp_sigma_mg_block

subroutine calculate_sc_chi_block(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    do j=blk_ylb,blk_yub
        do i=blk_xlb,blk_xub
            blk%chi_sc(i,j,1)=1d0
        end do
    end do
end subroutine calculate_sc_chi_block

subroutine calculate_rad_source_block(blk)
    type(blockdef), pointer :: blk
    real(8) :: temp
    integer :: i,j
    if (nd==1) then
    else if (nd==2) then
        do j=blk_ylb,blk_yub
            do i=blk_xlb,blk_xub
                temp=blk%temp(i,j,1)
                blk%rad_source(i,j,1)=blk%chi_sc(i,j,1)*planck_function(temp)
            end do
        end do
    end if
end subroutine calculate_rad_source_block

subroutine calculate_radpower_mg_block(blk)
    type(blockdef), pointer :: blk
    real(8) :: temp,radpower,dradpowerdT
    integer :: i,j
    if (nd==1) then
        do i=1,blk_size_nx
            do j=1,nmg
                temp=blk%temp(i,1,1)
                call radpower_mg(j,temp,radpower)
                call dradpowerdT_mg(j,temp,dradpowerdT)
                blk%a_mg(j,i,1,1)=radpower*4d0*pi/c_light
                blk%b_mg(j,i,1,1)=dradpowerdT*4d0*pi/c_light
            end do
        end do
    else
    end if
end subroutine calculate_radpower_mg_block

subroutine radpower_mg(igroup,temp,radpower)
    integer :: igroup,i
    real(8) :: temp,radpower,nu1,nu2,bnu1,bnu2,lnnu1,lnnu2
    radpower=0d0
    do i=1,nradcell
        nu1=rad_mg(i,igroup)
        nu2=rad_mg(i+1,igroup)
        lnnu1=log(nu1)
        lnnu2=log(nu2)
        bnu1=planck_law_frequency_dlnnu(temp,nu1)
        bnu2=planck_law_frequency_dlnnu(temp,nu2)
        radpower=radpower+(bnu1+bnu2)/2*(lnnu2-lnnu1)
    end do
end subroutine radpower_mg

subroutine dradpowerdT_mg(igroup,temp,dradpowerdT)
    integer :: igroup,i
    real(8) :: temp,dradpowerdT,nu1,nu2,dbnu1dT,dbnu2dT,lnnu1,lnnu2
    dradpowerdT=0d0
    do i=1,nradcell
        nu1=rad_mg(i,igroup)
        nu2=rad_mg(i+1,igroup)
        lnnu1=log(nu1)
        lnnu2=log(nu2)
        dbnu1dT=planck_law_dfrequency_dlnnu(temp,nu1)
        dbnu2dT=planck_law_dfrequency_dlnnu(temp,nu2)
        dradpowerdT=dradpowerdT+(dbnu1dT+dbnu2dT)/2*(lnnu2-lnnu1)
    end do
end subroutine dradpowerdT_mg

subroutine check_conductivity()
    type(blockdef), pointer :: blk
    if (rank==0) then
        blk=>blk_processor_tail
        print *,time_sys%ntimestep,blk%kx(blk_size_nx,1,1),blk%level
        print *,blk%Erad(blk_size_nx,1,1),blk%Erad_next,blk%sigma_rosseland(blk_size_nx,1,1),blk%sigma_rosseland_next
    end if
    if (rank==1) then
        blk=>blk_processor_head
        print *,time_sys%ntimestep,blk%kx(0,1,1),blk%level
        print *,blk%Erad_pre,blk%Erad(1,1,1),blk%sigma_rosseland_pre,blk%sigma_rosseland(1,1,1)
    end if
end subroutine check_conductivity

subroutine calculate_fld_conductivity()
    type(blockdef), pointer :: blk
    integer :: i
    if (nd==1) then
        call communicate_blocks_fld_1d()
        call calculate_fld_conductivity_1d()
    else if (nd==2) then
        call communicate_blocks_2d_fld()
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            select case (blk%loc_type)
            case (1)
                call calculate_fld_conductivity_2d_1(blk)
            case (2)
                call calculate_fld_conductivity_2d_2(blk)
            case (3)
                call calculate_fld_conductivity_2d_3(blk)
            case (4)
                call calculate_fld_conductivity_2d_4(blk)
            case (5)
                call calculate_fld_conductivity_2d_5(blk)
            case (6)
                call calculate_fld_conductivity_2d_6(blk)
            case (7)
                call calculate_fld_conductivity_2d_7(blk)
            case (8)
                call calculate_fld_conductivity_2d_8(blk)
            case (9)
                call calculate_fld_conductivity_2d_9(blk)
            end select
            if (i/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
        nullify(blk)
    end if
end subroutine calculate_fld_conductivity

subroutine calculate_fld_conductivity_1d()
    type(blockdef), pointer :: blk
    real(8) :: x1,x2,xi,dx,Erad_l,Erad_u,sigma_l,sigma_u
    integer :: i,j
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        dx=blk%dxyz(1)
        !first cell of this block
        if (associated(blk,blk_head)) then
            if (rad_bound_type(1)==2) then
                blk%kx(0,1,1)=c_light/blk%sigma_rosseland(1,1,1)/3d0
            else if (rad_bound_type(1)==9) then
            end if
        else
            x1=blk%nb_coor%xl(1)
            x2=blk%x_center(1)
            xi=blk%x_interface(0)
            Erad_l=blk%Erad_pre
            Erad_u=blk%Erad(1,1,1)
            sigma_l=blk%sigma_rosseland_pre
            sigma_u=blk%sigma_rosseland(1,1,1)
            blk%kx(0,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end if
        !interior cells
        do j=1,blk_size_nx-1
            x1=blk%x_center(j)
            x2=blk%x_center(j+1)
            xi=blk%x_interface(j)
            Erad_l=blk%Erad(j,1,1)
            Erad_u=blk%Erad(j+1,1,1)
            sigma_l=blk%sigma_rosseland(j,1,1)
            sigma_u=blk%sigma_rosseland(j+1,1,1)
            blk%kx(j,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
        !last cell of this block
        if (associated(blk,blk_tail)) then
            if (rad_bound_type(2)==2) then
                blk%kx(blk_size_nx,1,1)=c_light/blk%sigma_rosseland(blk_size_nx,1,1)/3d0
            else
                x1=blk%x_center(blk_size_nx)
                x2=blk%x_center(blk_size_nx+1)
                xi=blk%x_interface(blk_size_nx)
                Erad_l=blk%Erad(blk_size_nx,1,1)
                if (igeometry==0) then
                    Erad_u=blk%Erad(blk_size_nx,1,1)/(1d0+3d0*blk%sigma_rosseland(blk_size_nx,1,1)*dx/2d0)
                    blk%Erad(blk_size_nx+1,1,1)=Erad_u
                else if (igeometry==2) then
                    Erad_u=blk%Erad(blk_size_nx,1,1)*x1**2/x2**2
                    blk%Erad(blk_size_nx+1,1,1)=Erad_u
                end if
                sigma_l=blk%sigma_rosseland(blk_size_nx,1,1)
                sigma_u=blk%sigma_rosseland(blk_size_nx,1,1)
                blk%kx(blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
            end if
        else
            x1=blk%x_center(blk_size_nx)
            x2=blk%nb_coor%xu(1)
            xi=blk%x_interface(blk_size_nx)
            Erad_l=blk%Erad(blk_size_nx,1,1)
            Erad_u=blk%Erad_next
            sigma_l=blk%sigma_rosseland(blk_size_nx,1,1)
            sigma_u=blk%sigma_rosseland_next
            blk%kx(blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    nullify(blk)
end subroutine calculate_fld_conductivity_1d

!subroutine calculate_fld_conductivity_xl_1d_block(blk)
!    !calculate the xl block of a processor
!    type(blockdef), pointer :: blk
!    real(8) :: Erad_l,Erad_u,sigma_l,sigma_u,dx,dx_pre,dxyz_block(3),pos_block(3),Erad,sigma
!    integer :: i,ijk(3)
!    !calculate the conductivity of the xl interface
!    if (rank==0) then
!        if (rad_bound_type(1)==2) then
!            blk%kx(0,1,1)=c_light/blk%sigma_rosseland(1,1,1)/3d0
!        else if (rad_bound_type(1)==9) then
!            dx=blk%dxyz(1)
!            Erad_l=blk%Erad(0,1,1)
!            Erad_u=blk%Erad(1,1,1)
!            sigma_l=blk%sigma_rosseland(0,1,1)
!            sigma_u=blk%sigma_rosseland(1,1,1)
!            blk%kx(0,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!        end if
!    else
!        dx=(blk%dxyz(1)+blk%dx_pre)/2
!        Erad_l=blk%Erad_pre
!        Erad_u=blk%Erad(1,1,1)
!        sigma_l=blk%sigma_rosseland_pre
!        sigma_u=blk%sigma_rosseland(1,1,1)
!        blk%kx(0,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    end if
!    !calculate the conductivity of the interior interfaces
!    dx=blk%dxyz(1)
!    do i=1,blk_size_nx-1
!        Erad_l=blk%Erad(i,1,1)
!        Erad_u=blk%Erad(i+1,1,1)
!        sigma_l=blk%sigma_rosseland(i,1,1)
!        sigma_u=blk%sigma_rosseland(i+1,1,1)
!        blk%kx(i,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    end do
!    !calculate the conductivity of the xu interface
!    dx=(blk%dxyz(1)+blk%dx_next)/2
!    Erad_l=blk%Erad(blk_size_nx,1,1)
!    Erad_u=blk%Erad_next
!    sigma_l=blk%sigma_rosseland(blk_size_nx,1,1)
!    sigma_u=blk%sigma_rosseland_next
!    blk%kx(blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!end subroutine calculate_fld_conductivity_xl_1d_block
!
!subroutine calculate_fld_conductivity_xu_1d_block(blk)
!    !calculate the xu block of a processor
!    type(blockdef), pointer :: blk
!    real(8) :: Erad_l,Erad_u,sigma_l,sigma_u,dx,dx_next,r,r_out,rthetaphi(3),dxyz_block(3),pos_block(3),Erad,sigma
!    integer :: i,ijk(3)
!    !calculate the conductivity of the xl interface
!    dx=(blk%dxyz(1)+blk%dx_pre)/2
!    Erad_l=blk%Erad_pre
!    Erad_u=blk%Erad(1,1,1)
!    sigma_l=blk%sigma_rosseland_pre
!    sigma_u=blk%sigma_rosseland(1,1,1)
!    blk%kx(0,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    !calculate the conductivity of the interior interfaces
!    dx=blk%dxyz(1)
!    do i=1,blk_size_nx-1
!        Erad_l=blk%Erad(i,1,1)
!        Erad_u=blk%Erad(i+1,1,1)
!        sigma_l=blk%sigma_rosseland(i,1,1)
!        sigma_u=blk%sigma_rosseland(i+1,1,1)
!        blk%kx(i,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    end do
!    !calculate the conductivity of the xu interface
!    if (rank==np-1) then
!        if (rad_bound_type(2)==1) then
!            if (igeometry==0) then
!                Erad_l=blk%Erad(blk_size_nx,1,1)
!                Erad_u=blk%Erad(blk_size_nx,1,1)/(1d0+3d0*blk%sigma_rosseland(blk_size_nx,1,1)*dx/2d0)
!                blk%Erad(blk_size_nx+1,1,1)=Erad_u
!                sigma_l=blk%sigma_rosseland(blk_size_nx,1,1)
!                sigma_u=blk%sigma_rosseland(blk_size_nx,1,1)
!                blk%kx(blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!            else if (igeometry==2) then
!                dxyz_block=blk%dxyz
!                pos_block=blk%pos
!                ijk=(/blk_size_nx,1,1/)
!                call ijk_to_coords(ijk,dxyz_block,pos_block,rthetaphi)
!                r=rthetaphi(1)
!                r_out=r+dx
!                Erad_l=blk%Erad(blk_size_nx,1,1)
!                Erad_u=blk%Erad(blk_size_nx,1,1)*r**2/r_out**2
!                blk%Erad(blk_size_nx+1,1,1)=Erad_u
!                sigma_l=blk%sigma_rosseland(blk_size_nx,1,1)
!                sigma_u=blk%sigma_rosseland(blk_size_nx,1,1)
!                blk%kx(blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!            end if
!        else if (rad_bound_type(2)==2) then
!            blk%kx(blk_size_nx,1,1)=c_light/blk%sigma_rosseland(blk_size_nx,1,1)/3d0
!        else if (rad_bound_type(2)==9) then
!            Erad_l=blk%Erad(blk_size_nx,1,1)
!            Erad_u=blk%Erad(blk_size_nx+1,1,1)
!            sigma_l=blk%sigma_rosseland(blk_size_nx,1,1)
!            sigma_u=blk%sigma_rosseland(blk_size_nx+1,1,1)
!            blk%kx(blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!        end if
!    else
!        dx=(blk%dxyz(1)+blk%dx_next)/2
!        Erad_l=blk%Erad(blk_size_nx,1,1)
!        Erad_u=blk%Erad_next
!        sigma_l=blk%sigma_rosseland(blk_size_nx,1,1)
!        sigma_u=blk%sigma_rosseland_next
!        blk%kx(blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    end if
!end subroutine calculate_fld_conductivity_xu_1d_block
!
!subroutine calculate_fld_conductivity_interior_1d_block(blk)
!    !calculate the conductivity of blocks that are interior to a processor
!    type(blockdef), pointer :: blk
!    real(8) :: Erad_l,Erad_u,sigma_l,sigma_u,dx,dxyz_block(3),pos_block(3)
!    integer :: i,ijk(3)
!    !calculate the conductivity of the xl interface
!    dx=(blk%dxyz(1)+blk%dx_pre)/2
!    Erad_l=blk%Erad_pre
!    Erad_u=blk%Erad(1,1,1)
!    sigma_l=blk%sigma_rosseland_pre
!    sigma_u=blk%sigma_rosseland(1,1,1)
!    blk%kx(0,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    !calculate the conductivity of the interior interfaces
!    dx=blk%dxyz(1)
!    do i=1,blk_size_nx-1
!        Erad_l=blk%Erad(i,1,1)
!        Erad_u=blk%Erad(i+1,1,1)
!        sigma_l=blk%sigma_rosseland(i,1,1)
!        sigma_u=blk%sigma_rosseland(i+1,1,1)
!        blk%kx(i,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    end do
!    !calculate the conductivity of the xu interface
!    dx=(blk%dxyz(1)+blk%dx_next)/2
!    Erad_l=blk%Erad(blk_size_nx,1,1)
!    Erad_u=blk%Erad_next
!    sigma_l=blk%sigma_rosseland(blk_size_nx,1,1)
!    sigma_u=blk%sigma_rosseland_next
!    blk%kx(blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!end subroutine calculate_fld_conductivity_interior_1d_block

subroutine calculate_fld_conductivity_2d_1(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dx,Erad_l,Erad_u,Erad_ext,sigma_l,sigma_u,sigma_ext,x1,x2,xi,y1,y2,yi
    dx=blk%dxyz(1)
    !x direction
    if (rad_bound_type(1)==2) then
        blk%kx(0,1:blk_size_ny,1)=c_light/blk%sigma_rosseland(1,1:blk_size_ny,1)/3d0
    else if (rad_bound_type(1)==1) then
        do j=1,blk_size_ny
            Erad_u=blk%Erad(1,j,1)
            Erad_ext=blk%Erad(1,j,1)/(1d0+3d0*blk%sigma_rosseland(1,j,1)*dx/2d0)
            sigma_u=blk%sigma_rosseland(1,j,1)
            sigma_ext=sigma_u
            xi=blk%x_interface(0)
            x1=xi-dx/2d0
            x2=blk%x_center(1)
            blk%kx(0,j,1)=fld_conductivity(Erad_ext,Erad_u,sigma_ext,sigma_u,x1,x2,xi)
        end do
    else
    end if
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i+1,j,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i+1,j,1)
            xi=blk%x_interface(i)
            x1=blk%x_center(i)
            x2=blk%x_center(i+1)
            blk%kx(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
    end do
    do j=1,blk_size_ny
        Erad_l=blk%Erad(blk_size_nx,j,1)
        Erad_u=blk%Erad_xu(j)
        sigma_l=blk%sigma_rosseland(blk_size_nx,j,1)
        sigma_u=blk%sigma_rosseland_xu(j)
        xi=blk%x_interface(blk_size_nx)
        x1=blk%x_center(blk_size_nx)
        x2=xi+dx/2
        blk%kx(blk_size_nx,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    !y direction
    do i=1,blk_size_nx
        Erad_l=blk%Erad_yl(i)
        Erad_u=blk%Erad(i,1,1)
        sigma_l=blk%sigma_rosseland_yl(i)
        sigma_u=blk%sigma_rosseland(i,1,1)
        yi=blk%y_interface(0)
        y1=yi-dx/2
        y2=blk%y_center(1)
        blk%ky(i,0,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i,j+1,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i,j+1,1)
            yi=blk%y_interface(j)
            y1=blk%y_center(j)
            y2=blk%y_center(j+1)
            blk%ky(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
        end do
    end do
    if (rad_bound_type(4)==2) then
        blk%ky(1:blk_size_nx,blk_size_ny,1)=c_light/blk%sigma_rosseland(1:blk_size_nx,blk_size_ny,1)/3d0
    else if (rad_bound_type(4)==1) then
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,blk_size_ny,1)
            Erad_ext=blk%Erad(i,blk_size_ny,1)/(1d0+3d0*blk%sigma_rosseland(i,blk_size_ny,1)*dx/2d0)
            sigma_l=blk%sigma_rosseland(i,blk_size_ny,1)
            sigma_ext=sigma_l
            yi=blk%y_interface(blk_size_ny)
            y1=blk%y_center(blk_size_ny)
            y2=yi+dx/2
            blk%ky(i,blk_size_ny,1)=fld_conductivity(Erad_l,Erad_ext,sigma_l,sigma_ext,y1,y2,yi)
        end do
    else
    end if
end subroutine calculate_fld_conductivity_2d_1

subroutine calculate_fld_conductivity_2d_2(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dx,Erad_l,Erad_u,Erad_ext,sigma_l,sigma_u,sigma_ext,x1,x2,xi,y1,y2,yi
    dx=blk%dxyz(1)
    !x direction
    do j=1,blk_size_ny
        Erad_l=blk%Erad_xl(j)
        Erad_u=blk%Erad(1,j,1)
        sigma_l=blk%sigma_rosseland_xl(j)
        sigma_u=blk%sigma_rosseland(1,j,1)
        xi=blk%x_interface(0)
        x1=xi-dx/2
        x2=blk%x_center(1)
        blk%kx(0,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i+1,j,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i+1,j,1)
            xi=blk%x_interface(i)
            x1=blk%x_center(i)
            x2=blk%x_center(i+1)
            blk%kx(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
    end do
    do j=1,blk_size_ny
        Erad_l=blk%Erad(blk_size_nx,j,1)
        Erad_u=blk%Erad_xu(j)
        sigma_l=blk%sigma_rosseland(blk_size_nx,j,1)
        sigma_u=blk%sigma_rosseland_xu(j)
        xi=blk%x_interface(blk_size_nx)
        x1=blk%x_center(blk_size_nx)
        x2=xi+dx/2
        blk%kx(blk_size_nx,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    !y direction
    do i=1,blk_size_nx
        Erad_l=blk%Erad_yl(i)
        Erad_u=blk%Erad(i,1,1)
        sigma_l=blk%sigma_rosseland_yl(i)
        sigma_u=blk%sigma_rosseland(i,1,1)
        yi=blk%y_interface(0)
        y1=yi-dx/2
        y2=blk%y_center(1)
        blk%ky(i,0,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i,j+1,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i,j+1,1)
            yi=blk%y_interface(j)
            y1=blk%y_center(j)
            y2=blk%y_center(j+1)
            blk%ky(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
        end do
    end do
    if (rad_bound_type(4)==2) then
        blk%ky(1:blk_size_nx,blk_size_ny,1)=c_light/blk%sigma_rosseland(1:blk_size_nx,blk_size_ny,1)/3d0
    else if (rad_bound_type(4)==1) then
    else
    end if
end subroutine calculate_fld_conductivity_2d_2

subroutine calculate_fld_conductivity_2d_3(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dx,Erad_l,Erad_u,Erad_ext,sigma_l,sigma_u,sigma_ext,x1,x2,xi,y1,y2,yi
    dx=blk%dxyz(1)
    !x direction
    do j=1,blk_size_ny
        Erad_l=blk%Erad_xl(j)
        Erad_u=blk%Erad(1,j,1)
        sigma_l=blk%sigma_rosseland_xl(j)
        sigma_u=blk%sigma_rosseland(1,j,1)
        xi=blk%x_interface(0)
        x1=xi-dx/2
        x2=blk%x_center(1)
        blk%kx(0,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i+1,j,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i+1,j,1)
            xi=blk%x_interface(i)
            x1=blk%x_center(i)
            x2=blk%x_center(i+1)
            blk%kx(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
    end do
    if (rad_bound_type(2)==2) then
        blk%kx(blk_size_nx,1:blk_size_ny,1)=c_light/blk%sigma_rosseland(blk_size_nx,1:blk_size_ny,1)/3d0
    else if (rad_bound_type(2)==1) then
    else
    end if
    !y direction
    do i=1,blk_size_nx
        Erad_l=blk%Erad_yl(i)
        Erad_u=blk%Erad(i,1,1)
        sigma_l=blk%sigma_rosseland_yl(i)
        sigma_u=blk%sigma_rosseland(i,1,1)
        yi=blk%y_interface(0)
        y1=yi-dx/2
        y2=blk%y_center(1)
        blk%ky(i,0,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i,j+1,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i,j+1,1)
            yi=blk%y_interface(j)
            y1=blk%y_center(j)
            y2=blk%y_center(j+1)
            blk%ky(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
        end do
    end do
    if (rad_bound_type(4)==2) then
        blk%ky(1:blk_size_nx,blk_size_ny,1)=c_light/blk%sigma_rosseland(1:blk_size_nx,blk_size_ny,1)/3d0
    else if (rad_bound_type(4)==1) then
    else
    end if
end subroutine calculate_fld_conductivity_2d_3

subroutine calculate_fld_conductivity_2d_4(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dx,Erad_l,Erad_u,Erad_ext,sigma_l,sigma_u,sigma_ext,x1,x2,xi,y1,y2,yi
    dx=blk%dxyz(1)
    !x direction
    if (rad_bound_type(1)==2) then
        blk%kx(0,1:blk_size_ny,1)=c_light/blk%sigma_rosseland(0,1:blk_size_ny,1)/3d0
    else if (rad_bound_type(1)==1) then
    else
    end if
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i+1,j,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i+1,j,1)
            xi=blk%x_interface(i)
            x1=blk%x_center(i)
            x2=blk%x_center(i+1)
            blk%kx(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
    end do
    do j=1,blk_size_ny
        Erad_l=blk%Erad(blk_size_nx,j,1)
        Erad_u=blk%Erad_xu(j)
        sigma_l=blk%sigma_rosseland(blk_size_nx,j,1)
        sigma_u=blk%sigma_rosseland_xu(j)
        xi=blk%x_interface(blk_size_nx)
        x1=blk%x_center(blk_size_nx)
        x2=xi+dx/2
        blk%kx(blk_size_nx,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    !y direction
    do i=1,blk_size_nx
        Erad_l=blk%Erad_yl(i)
        Erad_u=blk%Erad(i,1,1)
        sigma_l=blk%sigma_rosseland_yl(i)
        sigma_u=blk%sigma_rosseland(i,1,1)
        yi=blk%y_interface(0)
        y1=yi-dx/2
        y2=blk%y_center(1)
        blk%ky(i,0,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i,j+1,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i,j+1,1)
            yi=blk%y_interface(j)
            y1=blk%y_center(j)
            y2=blk%y_center(j+1)
            blk%ky(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
        end do
    end do
    do i=1,blk_size_nx
        Erad_l=blk%Erad(i,blk_size_ny,1)
        Erad_u=blk%Erad_yu(i)
        sigma_l=blk%sigma_rosseland(i,blk_size_ny,1)
        sigma_u=blk%sigma_rosseland_yu(i)
        yi=blk%y_interface(blk_size_ny)
        y1=blk%y_center(blk_size_ny)
        y2=yi+dx/2
        blk%ky(i,blk_size_ny,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
end subroutine calculate_fld_conductivity_2d_4

subroutine calculate_fld_conductivity_2d_5(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dx,Erad_l,Erad_u,Erad_ext,sigma_l,sigma_u,sigma_ext,x1,x2,xi,y1,y2,yi
    dx=blk%dxyz(1)
    !x direction
    do j=1,blk_size_ny
        Erad_l=blk%Erad_xl(j)
        Erad_u=blk%Erad(1,j,1)
        sigma_l=blk%sigma_rosseland_xl(j)
        sigma_u=blk%sigma_rosseland(1,j,1)
        xi=blk%x_interface(0)
        x1=xi-dx/2d0
        x2=blk%x_center(1)
        blk%kx(0,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i+1,j,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i+1,j,1)
            xi=blk%x_interface(i)
            x1=blk%x_center(i)
            x2=blk%x_center(i+1)
            blk%kx(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
    end do
    do j=1,blk_size_ny
        Erad_l=blk%Erad(blk_size_nx,j,1)
        Erad_u=blk%Erad_xu(j)
        sigma_l=blk%sigma_rosseland(blk_size_nx,j,1)
        sigma_u=blk%sigma_rosseland_xu(j)
        xi=blk%x_interface(blk_size_nx)
        x1=blk%x_center(blk_size_nx)
        x2=xi+dx/2
        blk%kx(blk_size_nx,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    !y direction
    do i=1,blk_size_nx
        Erad_l=blk%Erad_yl(i)
        Erad_u=blk%Erad(i,1,1)
        sigma_l=blk%sigma_rosseland_yl(i)
        sigma_u=blk%sigma_rosseland(i,1,1)
        yi=blk%y_interface(0)
        y1=yi-dx/2
        y2=blk%y_center(1)
        blk%ky(i,0,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i,j+1,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i,j+1,1)
            yi=blk%y_interface(j)
            y1=blk%y_center(j)
            y2=blk%y_center(j+1)
            blk%ky(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
        end do
    end do
    do i=1,blk_size_nx
        Erad_l=blk%Erad(i,blk_size_ny,1)
        Erad_u=blk%Erad_yu(i)
        sigma_l=blk%sigma_rosseland(i,blk_size_ny,1)
        sigma_u=blk%sigma_rosseland_yu(i)
        yi=blk%y_interface(blk_size_ny)
        y1=blk%y_center(blk_size_ny)
        y2=yi+dx/2
        blk%ky(i,blk_size_ny,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
end subroutine calculate_fld_conductivity_2d_5

subroutine calculate_fld_conductivity_2d_6(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dx,Erad_l,Erad_u,Erad_ext,sigma_l,sigma_u,sigma_ext,x1,x2,xi,y1,y2,yi
    dx=blk%dxyz(1)
    !x direction
    do j=1,blk_size_ny
        Erad_l=blk%Erad_xl(j)
        Erad_u=blk%Erad(1,j,1)
        sigma_l=blk%sigma_rosseland_xl(j)
        sigma_u=blk%sigma_rosseland(1,j,1)
        xi=blk%x_interface(0)
        x1=xi-dx/2d0
        x2=blk%x_center(1)
        blk%kx(0,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i+1,j,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i+1,j,1)
            xi=blk%x_interface(i)
            x1=blk%x_center(i)
            x2=blk%x_center(i+1)
            blk%kx(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
    end do
    if (rad_bound_type(2)==2) then
        blk%kx(blk_size_nx,1:blk_size_ny,1)=c_light/blk%sigma_rosseland(blk_size_nx,1:blk_size_ny,1)/3d0
    else if (rad_bound_type(2)==1) then
    else
    end if
    !y direction
    do i=1,blk_size_nx
        Erad_l=blk%Erad_yl(i)
        Erad_u=blk%Erad(i,1,1)
        sigma_l=blk%sigma_rosseland_yl(i)
        sigma_u=blk%sigma_rosseland(i,1,1)
        yi=blk%y_interface(0)
        y1=yi-dx/2
        y2=blk%y_center(1)
        blk%ky(i,0,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i,j+1,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i,j+1,1)
            yi=blk%y_interface(j)
            y1=blk%y_center(j)
            y2=blk%y_center(j+1)
            blk%ky(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
        end do
    end do
    do i=1,blk_size_nx
        Erad_l=blk%Erad(i,blk_size_ny,1)
        Erad_u=blk%Erad_yu(i)
        sigma_l=blk%sigma_rosseland(i,blk_size_ny,1)
        sigma_u=blk%sigma_rosseland_yu(i)
        yi=blk%y_interface(blk_size_ny)
        y1=blk%y_center(blk_size_ny)
        y2=yi+dx/2
        blk%ky(i,blk_size_ny,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
end subroutine calculate_fld_conductivity_2d_6

subroutine calculate_fld_conductivity_2d_7(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dx,Erad_l,Erad_u,Erad_ext,sigma_l,sigma_u,sigma_ext,x1,x2,xi,y1,y2,yi
    dx=blk%dxyz(1)
    !x direction
    if (rad_bound_type(1)==2) then
        blk%kx(0,1:blk_size_ny,1)=c_light/blk%sigma_rosseland(1,1:blk_size_ny,1)/3d0
    else if (rad_bound_type(1)==1) then
    else
    end if
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i+1,j,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i+1,j,1)
            xi=blk%x_interface(i)
            x1=blk%x_center(i)
            x2=blk%x_center(i+1)
            blk%kx(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
    end do
    do j=1,blk_size_ny
        Erad_l=blk%Erad(blk_size_nx,j,1)
        Erad_u=blk%Erad_xu(j)
        sigma_l=blk%sigma_rosseland(blk_size_nx,j,1)
        sigma_u=blk%sigma_rosseland_xu(j)
        xi=blk%x_interface(blk_size_nx)
        x1=blk%x_center(blk_size_nx)
        x2=xi+dx/2
        blk%kx(blk_size_nx,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    !y direction
    if (rad_bound_type(3)==2) then
        blk%ky(1:blk_size_nx,0,1)=c_light/blk%sigma_rosseland(1:blk_size_nx,1,1)/3d0
    else if (rad_bound_type(3)==1) then
    else
    end if
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i,j+1,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i,j+1,1)
            yi=blk%y_interface(j)
            y1=blk%y_center(j)
            y2=blk%y_center(j+1)
            blk%ky(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
        end do
    end do
    do i=1,blk_size_nx
        Erad_l=blk%Erad(i,blk_size_ny,1)
        Erad_u=blk%Erad_yu(i)
        sigma_l=blk%sigma_rosseland(i,blk_size_ny,1)
        sigma_u=blk%sigma_rosseland_yu(i)
        yi=blk%y_interface(blk_size_ny)
        y1=blk%y_center(blk_size_ny)
        y2=yi+dx/2
        blk%ky(i,blk_size_ny,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
end subroutine calculate_fld_conductivity_2d_7

subroutine calculate_fld_conductivity_2d_8(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dx,Erad_l,Erad_u,Erad_ext,sigma_l,sigma_u,sigma_ext,x1,x2,xi,y1,y2,yi
    dx=blk%dxyz(1)
    !x direction
    do j=1,blk_size_ny
        Erad_l=blk%Erad_xl(j)
        Erad_u=blk%Erad(1,j,1)
        sigma_l=blk%sigma_rosseland_xl(j)
        sigma_u=blk%sigma_rosseland(1,j,1)
        xi=blk%x_interface(0)
        x1=xi-dx/2d0
        x2=blk%x_center(1)
        blk%kx(0,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i+1,j,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i+1,j,1)
            xi=blk%x_interface(i)
            x1=blk%x_center(i)
            x2=blk%x_center(i+1)
            blk%kx(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
    end do
    do j=1,blk_size_ny
        Erad_l=blk%Erad(blk_size_nx,j,1)
        Erad_u=blk%Erad_xu(j)
        sigma_l=blk%sigma_rosseland(blk_size_nx,j,1)
        sigma_u=blk%sigma_rosseland_xu(j)
        xi=blk%x_interface(blk_size_nx)
        x1=blk%x_center(blk_size_nx)
        x2=xi+dx/2
        blk%kx(blk_size_nx,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    !y direction
    if (rad_bound_type(3)==2) then
        blk%ky(1:blk_size_nx,0,1)=c_light/blk%sigma_rosseland(1:blk_size_nx,1,1)/3d0
    else if (rad_bound_type(3)==1) then
    else
    end if
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i,j+1,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i,j+1,1)
            yi=blk%y_interface(j)
            y1=blk%y_center(j)
            y2=blk%y_center(j+1)
            blk%ky(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
        end do
    end do
    do i=1,blk_size_nx
        Erad_l=blk%Erad(i,blk_size_ny,1)
        Erad_u=blk%Erad_yu(i)
        sigma_l=blk%sigma_rosseland(i,blk_size_ny,1)
        sigma_u=blk%sigma_rosseland_yu(i)
        yi=blk%y_interface(blk_size_ny)
        y1=blk%y_center(blk_size_ny)
        y2=yi+dx/2
        blk%ky(i,blk_size_ny,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
end subroutine calculate_fld_conductivity_2d_8

subroutine calculate_fld_conductivity_2d_9(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: dx,Erad_l,Erad_u,Erad_ext,sigma_l,sigma_u,sigma_ext,x1,x2,xi,y1,y2,yi
    dx=blk%dxyz(1)
    !x direction
    do j=1,blk_size_ny
        Erad_l=blk%Erad_xl(j)
        Erad_u=blk%Erad(1,j,1)
        sigma_l=blk%sigma_rosseland_xl(j)
        sigma_u=blk%sigma_rosseland(1,j,1)
        xi=blk%x_interface(0)
        x1=xi-dx/2d0
        x2=blk%x_center(1)
        blk%kx(0,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    end do
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i+1,j,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i+1,j,1)
            xi=blk%x_interface(i)
            x1=blk%x_center(i)
            x2=blk%x_center(i+1)
            blk%kx(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
        end do
    end do
    if (rad_bound_type(2)==2) then
        blk%kx(blk_size_nx,1:blk_size_ny,1)=c_light/blk%sigma_rosseland(blk_size_nx,1:blk_size_ny,1)/3d0
    else if (rad_bound_type(2)==1) then
    else
    end if
    !y direction
    if (rad_bound_type(3)==2) then
        blk%ky(1:blk_size_nx,0,1)=c_light/blk%sigma_rosseland(1:blk_size_nx,1,1)/3d0
    else if (rad_bound_type(3)==1) then
    else
    end if
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            Erad_l=blk%Erad(i,j,1)
            Erad_u=blk%Erad(i,j+1,1)
            sigma_l=blk%sigma_rosseland(i,j,1)
            sigma_u=blk%sigma_rosseland(i,j+1,1)
            yi=blk%y_interface(j)
            y1=blk%y_center(j)
            y2=blk%y_center(j+1)
            blk%ky(i,j,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
        end do
    end do
    do i=1,blk_size_nx
        Erad_l=blk%Erad(i,blk_size_ny,1)
        Erad_u=blk%Erad_yu(i)
        sigma_l=blk%sigma_rosseland(i,blk_size_ny,1)
        sigma_u=blk%sigma_rosseland_yu(i)
        yi=blk%y_interface(blk_size_ny)
        y1=blk%y_center(blk_size_ny)
        y2=yi+dx/2
        blk%ky(i,blk_size_ny,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,y1,y2,yi)
    end do
end subroutine calculate_fld_conductivity_2d_9

subroutine calculate_fld_conductivity_mg()
!    type(blockdef), pointer :: blk
!    integer :: i
!    call communicate_blocks_1d_fld_mg()
!    blk=>blk_processor_head
!    call calculate_fld_conductivity_mg_xl_1d_block(blk)
!    blk=>blk%blk_next
!    do i=2,np_nblk(rank+1)-1
!        call calculate_fld_conductivity_mg_interior_1d_block(blk)
!        blk=>blk%blk_next
!    end do
!    call calculate_fld_conductivity_mg_xu_1d_block(blk)
end subroutine calculate_fld_conductivity_mg
!
!subroutine calculate_fld_conductivity_mg_xl_1d_block(blk)
!    !calculate the xl block of a processor, multigroup version
!    type(blockdef), pointer :: blk
!    real(8) :: Erad_l,Erad_u,sigma_l,sigma_u,dx,dx_pre,dxyz_block(3),pos_block(3),Erad,sigma
!    integer :: i,j,ijk(3)
!    !calculate the conductivity of the xl interface
!    if (rank==0) then
!        if (rad_bound_type(1)==2) then
!            do j=1,nmg
!                blk%kx_mg(j,0,1,1)=c_light/blk%sigma_rosseland_mg(j,1,1,1)/3d0
!            end do
!        else if (rad_bound_type(1)==9) then
!            do j=1,nmg
!                dx=blk%dxyz(1)
!                Erad_l=blk%Erad_mg(j,0,1,1)
!                Erad_u=blk%Erad_mg(j,1,1,1)
!                sigma_l=blk%sigma_rosseland_mg(j,0,1,1)
!                sigma_u=blk%sigma_rosseland_mg(j,1,1,1)
!                blk%kx_mg(j,0,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!            end do
!        end if
!    else
!        do j=1,nmg
!            dx=(blk%dxyz(1)+blk%dx_pre)/2
!            Erad_l=blk%Erad_pre_mg(j)
!            Erad_u=blk%Erad_mg(j,1,1,1)
!            sigma_l=blk%sigma_rosseland_pre_mg(j)
!            sigma_u=blk%sigma_rosseland_mg(j,1,1,1)
!            blk%kx_mg(j,0,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!        end do
!    end if
!    !calcualte the conductivity of the interior interfaces
!    dx=blk%dxyz(1)
!    do i=1,blk_size_nx-1
!        do j=1,nmg
!            Erad_l=blk%Erad_mg(j,i,1,1)
!            Erad_u=blk%Erad_mg(j,i+1,1,1)
!            sigma_l=blk%sigma_rosseland_mg(j,i,1,1)
!            sigma_u=blk%sigma_rosseland_mg(j,i+1,1,1)
!            blk%kx_mg(j,i,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!        end do
!    end do
!    !calculate the conductivity of the xu interface
!    dx=(blk%dxyz(1)+blk%dx_next)/2
!    do j=1,nmg
!        Erad_l=blk%Erad_mg(j,blk_size_nx,1,1)
!        Erad_u=blk%Erad_next_mg(j)
!        sigma_l=blk%sigma_rosseland_mg(j,blk_size_nx,1,1)
!        sigma_u=blk%sigma_rosseland_next_mg(j)
!        blk%kx_mg(j,blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    end do
!end subroutine calculate_fld_conductivity_mg_xl_1d_block
!
!subroutine calculate_fld_conductivity_mg_xu_1d_block(blk)
!    !calculate the xu block of a processor
!    type(blockdef), pointer :: blk
!    real(8) :: Erad_l,Erad_u,sigma_l,sigma_u,dx,dx_next,r,r_out,rthetaphi(3),dxyz_block(3),pos_block(3),Erad,sigma
!    integer :: i,j,ijk(3)
!    !calculate the conductivity of the xl interface
!    dx=(blk%dxyz(1)+blk%dx_pre)/2
!    do j=1,nmg
!        Erad_l=blk%Erad_pre_mg(j)
!        Erad_u=blk%Erad_mg(j,1,1,1)
!        sigma_l=blk%sigma_rosseland_pre_mg(j)
!        sigma_u=blk%sigma_rosseland_mg(j,1,1,1)
!        blk%kx_mg(j,0,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    end do
!    !calculate the conductivity of the interior interfaces
!    dx=blk%dxyz(1)
!    do i=1,blk_size_nx-1
!        do j=1,nmg
!            Erad_l=blk%Erad_mg(j,i,1,1)
!            Erad_u=blk%Erad_mg(j,i+1,1,1)
!            sigma_l=blk%sigma_rosseland_mg(j,i,1,1)
!            sigma_u=blk%sigma_rosseland_mg(j,i+1,1,1)
!            blk%kx_mg(j,i,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!        end do
!    end do
!    !calculate the conductivity of the xu interface
!    if (rank==np-1) then
!        if (rad_bound_type(2)==1) then
!            if (igeometry==0) then
!                do j=1,nmg
!                    Erad_l=blk%Erad_mg(j,blk_size_nx,1,1)
!                    Erad_u=blk%Erad_mg(j,blk_size_nx,1,1)/(1d0+3d0*blk%sigma_rosseland_mg(j,blk_size_nx,1,1)*dx/2d0)
!                    blk%Erad_mg(j,blk_size_nx+1,1,1)=Erad_u
!                    sigma_l=blk%sigma_rosseland_mg(j,blk_size_nx,1,1)
!                    sigma_u=blk%sigma_rosseland_mg(j,blk_size_nx,1,1)
!                    blk%kx_mg(j,blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!                end do
!            else if (igeometry==2) then
!                dxyz_block=blk%dxyz
!                pos_block=blk%pos
!                ijk=(/blk_size_nx,1,1/)
!                call ijk_to_coords(ijk,dxyz_block,pos_block,rthetaphi)
!                r=rthetaphi(1)
!                r_out=r+dx
!                do j=1,nmg
!                    Erad_l=blk%Erad_mg(j,blk_size_nx,1,1)
!                    Erad_u=blk%Erad_mg(j,blk_size_nx,1,1)*r**2/r_out**2
!                    blk%Erad_mg(j,blk_size_nx+1,1,1)=Erad_u
!                    sigma_l=blk%sigma_rosseland_mg(j,blk_size_nx,1,1)
!                    sigma_u=blk%sigma_rosseland_mg(j,blk_size_nx,1,1)
!                    blk%kx_mg(j,blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!                end do
!            end if
!        else if (rad_bound_type(2)==2) then
!            do j=1,nmg
!                blk%kx_mg(j,blk_size_nx,1,1)=c_light/blk%sigma_rosseland_mg(j,blk_size_nx,1,1)/3d0
!            end do
!        else if (rad_bound_type(2)==9) then
!            do j=1,nmg
!                Erad_l=blk%Erad_mg(j,blk_size_nx,1,1)
!                Erad_u=blk%Erad_mg(j,blk_size_nx+1,1,1)
!                sigma_l=blk%sigma_rosseland_mg(j,blk_size_nx,1,1)
!                sigma_u=blk%sigma_rosseland_mg(j,blk_size_nx+1,1,1)
!                blk%kx_mg(j,blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!            end do
!        end if
!    else
!        dx=(blk%dxyz(1)+blk%dx_next)/2
!        do j=1,nmg
!            Erad_l=blk%Erad_mg(j,blk_size_nx,1,1)
!            Erad_u=blk%Erad_next_mg(j)
!            sigma_l=blk%sigma_rosseland_mg(j,blk_size_nx,1,1)
!            sigma_u=blk%sigma_rosseland_next_mg(j)
!            blk%kx_mg(j,blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!        end do
!    end if
!end subroutine calculate_fld_conductivity_mg_xu_1d_block
!
!subroutine calculate_fld_conductivity_mg_interior_1d_block(blk)
!    !calcualte the conductivity of blocks that are interior to a processor
!    type(blockdef), pointer :: blk
!    real(8) :: Erad_l,Erad_u,sigma_l,sigma_u,dx,r,r_out,rthetaphi(3),dxyz_block(3),pos_block(3)
!    integer :: i,j,ijk(3)
!    !calculate the conductivity of the xl interface
!    dx=(blk%dxyz(1)+blk%dx_pre)/2
!    do j=1,nmg
!        Erad_l=blk%Erad_pre_mg(j)
!        Erad_u=blk%Erad_mg(j,1,1,1)
!        sigma_l=blk%sigma_rosseland_pre_mg(j)
!        sigma_u=blk%sigma_rosseland_mg(j,1,1,1)
!        blk%kx_mg(j,0,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    end do
!    !calculate the conductivity of the interior interfaces
!    dx=blk%dxyz(1)
!    do i=1,blk_size_nx-1
!        do j=1,nmg
!            Erad_l=blk%Erad_mg(j,i,1,1)
!            Erad_u=blk%Erad_mg(j,i+1,1,1)
!            sigma_l=blk%sigma_rosseland_mg(j,i,1,1)
!            sigma_u=blk%sigma_rosseland_mg(j,i+1,1,1)
!            blk%kx_mg(j,i,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!        end do
!    end do
!    !calculate the conductivity of the xu interface
!    dx=(blk%dxyz(1)+blk%dx_next)/2
!    do j=1,nmg
!        Erad_l=blk%Erad_mg(j,blk_size_nx,1,1)
!        Erad_u=blk%Erad_next_mg(j)
!        sigma_l=blk%sigma_rosseland_mg(j,blk_size_nx,1,1)
!        sigma_u=blk%sigma_rosseland_next_mg(j)
!        blk%kx_mg(j,blk_size_nx,1,1)=fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,dx)
!    end do
!end subroutine calculate_fld_conductivity_mg_interior_1d_block

function fld_conductivity(Erad_l,Erad_u,sigma_l,sigma_u,x1,x2,xi)
    real(8) :: Erad_l,Erad_u,sigma_l,sigma_u,Erad_interface,sigma_interface,fld_conductivity
    real(8) :: Erad_grad,r,x1,x2,xi,dx1,dx2,dx,x
    dx1=xi-x1
    dx2=x2-xi
    dx=x2-x1
    x=dx1/dx
    sigma_interface=sigma_l+(sigma_u-sigma_l)*x
    Erad_interface=Erad_l+(Erad_u-Erad_l)*x
    Erad_grad=(Erad_u-Erad_l)/dx
    r=abs(Erad_grad)/sigma_interface/Erad_interface
    fld_conductivity=c_light*fld_flux_limiter(r)/sigma_interface
end function fld_conductivity

end module radiation_common_functions
