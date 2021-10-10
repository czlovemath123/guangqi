module eos_tabulated
use datastructure
use mathlib
use phylib
implicit none

!**************interp method, EOS Riemann solver related tables and quantities***********
type(table2d), public :: cs_trho,eg_trho,t_prho,p_trho,gamma_trho,eg_rhop,rho_egp,maw_trho,  &
    t_egrho,meta_arhop,meta_egmrhop,meta_gammarhop,meta_prhoegm
type(table1d), public :: gamma_t,egm_t,g_t,t_g,cs_t,maw_t,g_egm,egm_g,t_egm,cs_g
real(8), public :: rho_min_bound,rho_max_bound,p_min_bound,p_max_bound,cs_min_bound,cs_max_bound,  &
    t_min_bound,t_max_bound,gamma_max,gamma_min,g_min_bound,g_max_bound
!**************interp method, EOS Riemann solver related tables and quantities***********
contains

!calculate x direction flux in EOS condition from w
subroutine eos_wtoxflux(w,xflux)
    real(8), dimension(5) :: w
    real(8), dimension(5) :: xflux
    real(8) :: rho,vx,vy,vz,v,eg,p
    !eg is energy per unit mass
    rho=w(1)
    vx=w(2)
    vy=w(3)
    vz=w(4)
    v=sqrt(vx*vx+vy*vy+vz*vz)
    p=w(5)
    eg=eos_tabulated_internal_e_per_mass(rho,p)
    xflux(1)=rho*vx
    xflux(2)=rho*vx*vx+p
    xflux(3)=rho*vx*vy
    xflux(4)=rho*vx*vz
    xflux(5)=vx*(0.5*rho*v*v+eg*rho+p)
end subroutine eos_wtoxflux

subroutine eos_wtou(w,u)
    real(8), dimension(5) :: w
    real(8), dimension(5) :: u
    real(8) :: rho,vx,vy,vz,v,E,eg,p
    !for w, it is rho, vx, vy, vz, p
    !eg is energy per unit mass
    real(8) :: log10_eg,log10_rho,log10_temp,log10_p
    rho=w(1)
    vx=w(2)
    vy=w(3)
    vz=w(4)
    v=sqrt(vx*vx+vy*vy+vz*vz)
    p=w(5)
    eg=eos_tabulated_internal_e_per_mass(rho,p)
    E=rho*eg+0.5*rho*v*v
    u(1)=rho
    u(2)=rho*vx
    u(3)=rho*vy
    u(4)=rho*vz
    u(5)=E
end subroutine eos_wtou

subroutine eos_utow(u,w)
    real(8), dimension(5) :: w,u
    real(8) :: rho,vx,vy,vz,v,E,eg,egv,egm,p,p_log10,t,t_log10
    !for w, it is rho, vx, vy, vz, p
    !eg is energy per unit mass
    !egv is energy per unit volume
    rho=u(1)
    vx=u(2)/rho
    vy=u(3)/rho
    vz=u(4)/rho
    v=sqrt(vx*vx+vy*vy+vz*vz)
    egv=(u(5)-0.5*rho*v*v)
    p=eos_tabulated_pressure(rho,egv)
    w(1)=rho
    w(2)=vx
    w(3)=vy
    w(4)=vz
    w(5)=p
end subroutine eos_utow

subroutine rho_p_cs_boundary(m,n,o)
    integer :: m,n,o
    rho_min_bound=10**p_trho%ylist(1)
    rho_max_bound=10**p_trho%ylist(m)
    p_min_bound=t_prho%xlist(1)
    p_max_bound=t_prho%xlist(o)
    cs_min_bound=10**minval(cs_trho%t2d)
    cs_max_bound=10**maxval(cs_trho%t2d)
    t_min_bound=p_trho%xlist(1)
    t_max_bound=p_trho%xlist(n)
end subroutine rho_p_cs_boundary

subroutine generate_eos_tables(subdir1,subdir2)
    !generate all the EOS tables to be use
    !m,n,o are density, temp and pressure sampling points, respectively
    !cs,eg,p,t will all in log10log10log10 form
    !gamma is in log10log10linear form
    !if EOS_test is .true., the code is testing constant gamma RS
    integer :: m,n,o,i,j,k,nrho,np,negm
    character(len=30) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,subdir1,subdir2
    character(len=4) :: mark1,mark2,mark3,mark4,mark5,mark6
    real(8), dimension(:), allocatable :: temp1,temp2,rholist,tlist,plist,eglist,meta_rholist,meta_plist,  &
        meta_egmlist,arhop_data,egmrhop_data,gammarhop_data,prhoegm_data
    real(8), dimension(:,:), allocatable :: arhop_table,egmrhop_table,gammarhop_table,prhoegm_table,temp1_table
    real(8) :: pressure,gamm,rho,t_temp,rho_temp,p_temp,log10_rho,log10_temp,log10_p,  &
        bound(2),dlogrho,prho
    type(table2d) :: eg_trho_temp
    m=environment%table_dim(1)
    n=environment%table_dim(2)
    nrho=environment%meta_table_dim(1)
    np=environment%meta_table_dim(2)
    negm=environment%meta_table_dim(3)
    o=m+n
    bound=(/real(-10,kind=8),real(10,kind=8)/)
    write(mark1,'(I4)')m
    write(mark2,'(I4)')n
    write(mark3,'(I4)')o
    write(mark4,'(I4)')nrho
    write(mark5,'(I4)')np
    write(mark6,'(I4)')negm
    !for tabulated EoS gas
    fmt1=trim('(')//mark1//trim('ES18.8E2)')
    fmt2=trim('(')//mark2//trim('ES18.8E2)')
    fmt3=trim('(')//mark3//trim('ES18.8E2)')
    fmt4=trim('(')//mark4//trim('ES18.8E2)')
    fmt5=trim('(')//mark5//trim('ES18.8E2)')
    fmt6=trim('(')//mark6//trim('ES18.8E2)')
    open(unit=11,file='tables/realistic_gas/'//trim(subdir1)//trim(subdir2)//'/p.txt',status='old',action='read')
    open(unit=12,file='tables/realistic_gas/'//trim(subdir1)//trim(subdir2)//'/gamma.txt',status='old',action='read')
    open(unit=13,file='tables/realistic_gas/'//trim(subdir1)//trim(subdir2)//'/eg.txt',status='old',action='read')
    open(unit=14,file='tables/realistic_gas/'//trim(subdir1)//trim(subdir2)//'/maw.txt',status='old',action='read')
    open(unit=15,file='tables/realistic_gas/'//trim(subdir1)//trim(subdir2)//'/arhop.txt',status='old',action='read')
    open(unit=16,file='tables/realistic_gas/'//trim(subdir1)//trim(subdir2)//'/egmrhop.txt',status='old',action='read')
    open(unit=17,file='tables/realistic_gas/'//trim(subdir1)//trim(subdir2)//'/gammarhop.txt',status='old',action='read')
    open(unit=18,file='tables/realistic_gas/'//trim(subdir1)//trim(subdir2)//'/prhoegm.txt',status='old',action='read')
    call create_table_2d(n,m,p_trho)
    call create_table_2d(n,m,gamma_trho)
    call create_table_2d(n,m,cs_trho)
    call create_table_2d(n,m,eg_trho)
    call create_table_2d(n,m,eg_trho_temp)
    call create_table_2d(o,m,t_prho)
    call create_table_2d(m,n,eg_rhop)
    call create_table_2d(o,n,rho_egp)
    call create_table_2d(n,m,maw_trho)
    call create_table_2d(o,m,t_egrho)
    call create_table_2d(nrho,np,meta_arhop)
    call create_table_2d(nrho,np,meta_egmrhop)
    call create_table_2d(nrho,np,meta_gammarhop)
    call create_table_2d(nrho,negm,meta_prhoegm)
    allocate(rholist(m),tlist(n),plist(o),eglist(o),temp1(n*m),temp2(o*m),arhop_data(nrho*np),egmrhop_data(nrho*np),  &
        gammarhop_data(nrho*np),prhoegm_data(nrho*negm))
    allocate(meta_rholist(nrho),meta_plist(np),meta_egmlist(negm))
    allocate(arhop_table(np,nrho),egmrhop_table(np,nrho),gammarhop_table(np,nrho),prhoegm_table(negm,nrho),temp1_table(m,n))
    read(11,fmt1)rholist
    read(11,fmt2)tlist
    read(11,*)temp1
    do i=1,n
        p_trho%t2d(i,:)=temp1(1+(i-1)*m:i*m)
    end do
    p_trho%xlist=tlist
    p_trho%ylist=rholist
    call loglogloginversefunction(tlist,rholist,plist,p_trho%t2d,t_prho%t2d)
    t_prho%xlist=plist
    t_prho%ylist=rholist
    read(12,fmt1)rholist
    read(12,fmt2)tlist
    read(12,*)temp1
    do i=1,n
        gamma_trho%t2d(i,:)=temp1(1+(i-1)*m:i*m)
    end do
    gamma_trho%xlist=tlist
    gamma_trho%ylist=rholist
    read(13,fmt1)rholist
    read(13,fmt2)tlist
    read(13,*)temp1
    do i=1,n
        eg_trho%t2d(i,:)=temp1(1+(i-1)*m:i*m)
    end do
    eg_trho%xlist=tlist
    eg_trho%ylist=rholist
    eg_trho_temp=eg_trho
    call loglogloginversefunction(tlist,rholist,eglist,eg_trho_temp%t2d,t_egrho%t2d)
    t_egrho%xlist=eglist
    t_egrho%ylist=rholist
    read(14,fmt1)rholist
    read(14,fmt2)tlist
    read(14,*)temp1
    temp1_table=reshape(temp1,(/m,n/))
    do i=1,n
        maw_trho%t2d(i,:)=temp1(1+(i-1)*m:i*m)
    end do
    !print *,transpose(temp1_table)-maw_trho%t2d
    maw_trho%xlist=tlist
    maw_trho%ylist=rholist
    read(15,fmt4)meta_rholist
    read(15,fmt5)meta_plist
    read(15,*)arhop_data
    arhop_table=reshape(arhop_data,(/np,nrho/))
    meta_arhop%t2d=transpose(arhop_table)
    meta_arhop%xlist=meta_rholist
    meta_arhop%ylist=meta_plist
    read(16,fmt4)meta_rholist
    read(16,fmt5)meta_plist
    read(16,*)egmrhop_data
    egmrhop_table=reshape(egmrhop_data,(/np,nrho/))
    meta_egmrhop%t2d=transpose(egmrhop_table)
    meta_egmrhop%xlist=meta_rholist
    meta_egmrhop%ylist=meta_plist
    read(17,fmt4)meta_rholist
    read(17,fmt5)meta_plist
    read(17,*)gammarhop_data
    gammarhop_table=reshape(gammarhop_data,(/np,nrho/))
    meta_gammarhop%t2d=transpose(gammarhop_table)
    meta_gammarhop%xlist=meta_rholist
    meta_gammarhop%ylist=meta_plist
    read(18,fmt4)meta_rholist
    read(18,fmt6)meta_egmlist
    read(18,*)prhoegm_data
    prhoegm_table=reshape(prhoegm_data,(/negm,nrho/))
    meta_prhoegm%t2d=transpose(prhoegm_table)
    meta_prhoegm%xlist=meta_rholist
    meta_prhoegm%ylist=meta_egmlist
    do i=1,n
        do j=1,m
            pressure=p_trho%t2d(i,j)
            gamm=gamma_trho%t2d(i,j)
            rho=10**rholist(j)
            cs_trho%t2d(i,j)=log10(sqrt(gamm*pressure/rho))
        end do
    end do
    cs_trho%xlist=tlist
    cs_trho%ylist=rholist
    eg_trho%t2d=log10(eg_trho%t2d)
    p_trho%t2d=log10(p_trho%t2d)
    deallocate(tlist,rholist,plist,eglist,temp1,temp2)
    deallocate(meta_rholist,meta_plist,meta_egmlist,arhop_data,arhop_table,egmrhop_data,egmrhop_table,  &
        gammarhop_data,gammarhop_table,prhoegm_data,prhoegm_table,temp1_table)
    call rho_p_cs_boundary(m,n,o)
    do i=11,14
        close(i)
    end do
end subroutine generate_eos_tables

function eos_tabulated_sound_speed(rho,pressure)
    !for real gas, gamm=gamm(rho,p) is tabulated
    !the table should be 'eos_sound_speed.txt' in problem directory
    !cs is cs_trho table, t is t_prho table
    real(8) :: rho,pressure,t,log10_rho,log10_pressure,log10_cs,eos_tabulated_sound_speed,log10_temp
    logical :: radiation
    radiation=.false.
    log10_rho=log10(rho)
    log10_pressure=log10(pressure)
    if (isnan(log10_rho)) then
        print *,'cs, rho=',rho,pressure
    end if
    !call interpolation2d_linear(log10_pressure,log10_rho,log10_temp,t_prho%xlist,t_prho%ylist,t_prho%t2d)
    !call interpolation2d_linear(log10_temp,log10_rho,log10_cs,cs_trho%xlist,cs_trho%ylist,cs_trho%t2d)
    !eos_tabulated_sound_speed=10**log10_cs
    call interpolation2d_linear(log10_rho,log10_pressure,eos_tabulated_sound_speed,meta_arhop%xlist,meta_arhop%ylist,meta_arhop%t2d)
end function eos_tabulated_sound_speed

function eos_tabulated_internal_e_per_mass(rho,pressure)
    !for real gas, e=e(rho,p) is tabulated
    !the table should be 'eos_e.txt' in problem directory
    !e is eg_trho table, t is t_prho table
    real(8) :: rho,pressure,t,log10_rho,log10_pressure,log10_eg,log10_egm,eos_tabulated_internal_e_per_mass,log10_temp
    logical :: radiation
    radiation=.false.
    log10_rho=log10(rho)
    log10_pressure=log10(pressure)
    if (isnan(log10_rho)) then
        print *,'e, rho=',rho,pressure
    end if
    !call interpolation2d_linear(log10_pressure,log10_rho,log10_temp,t_prho%xlist,t_prho%ylist,t_prho%t2d)
    !call interpolation2d_linear(log10_temp,log10_rho,log10_eg,eg_trho%xlist,eg_trho%ylist,eg_trho%t2d)
    !eos_tabulated_internal_e_per_mass=(10**log10_eg)/rho
    call interpolation2d_linear(log10_rho,log10_pressure,log10_egm,meta_egmrhop%xlist,meta_egmrhop%ylist,meta_egmrhop%t2d)
    eos_tabulated_internal_e_per_mass=10**log10_egm
end function eos_tabulated_internal_e_per_mass

function eos_tabulated_pressure(rho,eg)
    !for real gas, p=p(rho,eg) is tabulated
    !eg is in the unit of erg per volume
    real(8) :: rho,eg,egv,t,eos_tabulated_pressure,log10_rho,log10_eg,log10_egm,log10_temp,log10_p
    logical :: radiation
    radiation=.false.
    log10_rho=log10(rho)
    log10_eg=log10(eg)
    log10_egm=log10(eg/rho)
    !call interpolation2d_linear(log10_eg,log10_rho,log10_temp,t_egrho%xlist,t_egrho%ylist,t_egrho%t2d)
    !call interpolation2d_linear(log10_temp,log10_rho,log10_p,p_trho%xlist,p_trho%ylist,p_trho%t2d)
    !eos_tabulated_pressure=10**log10_p
    call interpolation2d_linear(log10_rho,log10_egm,log10_p,meta_prhoegm%xlist,meta_prhoegm%ylist,meta_prhoegm%t2d)
    eos_tabulated_pressure=10**log10_p
end function eos_tabulated_pressure

function eos_tabulated_temp(rho,p)
    !for real gas, t=t(p,rho) is tabulated
    !t is in unit of K
    real(8) :: rho,p,eos_tabulated_temp,log10_rho,log10_p,log10_temp
    logical :: radiation
    radiation=.false.
    log10_rho=log10(rho)
    log10_p=log10(p)
    call interpolation2d_linear(log10_p,log10_rho,log10_temp,t_prho%xlist,t_prho%ylist,t_prho%t2d)
    eos_tabulated_temp=10**log10_temp
end function eos_tabulated_temp

function eos_tabulated_gamma(rho,p)
    !given rho and p, calculate Gamma. used in eos_hllc_tabulated
    real(8) :: rho,p,eos_tabulated_gamma,log10_rho,log10_p
    log10_rho=log10(rho)
    log10_p=log10(p)
    call interpolation2d_linear(log10_rho,log10_p,eos_tabulated_gamma,meta_gammarhop%xlist,meta_gammarhop%ylist,meta_gammarhop%t2d)
end function eos_tabulated_gamma

function eos_tabulated_egv(rho,temp)
    real(8) :: rho,temp,log10_rho,log10_eg,eos_tabulated_egv,log10_temp
    log10_rho=log10(rho)
    log10_temp=log10(temp)
    if (isnan(log10_temp)) then
        print *,'rho, temp=',rho,temp
    end if
    call interpolation2d_linear(log10_temp,log10_rho,log10_eg,eg_trho%xlist,eg_trho%ylist,eg_trho%t2d)
    eos_tabulated_egv=10**log10_eg
end function eos_tabulated_egv

function eos_tabulated_maw_2d(rho,temp)
    !calculate the mean atomic weight
    real(8) :: rho,temp,eos_tabulated_maw_2d,log10_rho,log10_temp
    log10_rho=log10(rho)
    log10_temp=log10(temp)
    if (isnan(log10_temp)) then
        print *,'rho, temp=',rho,temp
    end if
    call interpolation2d_linear(log10_temp,log10_rho,eos_tabulated_maw_2d,maw_trho%xlist,maw_trho%ylist,maw_trho%t2d)
end function eos_tabulated_maw_2d

end module eos_tabulated
