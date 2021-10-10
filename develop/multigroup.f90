program multigroup
implicit none
real(8), parameter  ::    &
    gr=6.67259d-8,              &   !gravitational constant
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
real(8), dimension(:), allocatable :: group_frequency,radpower,dradpowerdT
real(8), dimension(:,:), allocatable :: rad_mg
real(8) :: tmax,tmin,numax,numin,dlnnu_group,temp,pi
integer :: i,groupnumber,nradcell
pi=3.1415926
tmax=1d5
tmin=1d2
groupnumber=10
nradcell=20
allocate(group_frequency(groupnumber+1),radpower(groupnumber),dradpowerdT(groupnumber))
numax=planck_function_peak_frequency(tmax)
numin=planck_function_peak_frequency(tmin)
dlnnu_group=log(numax/numin)/groupnumber
group_frequency(1)=numin
do i=2,groupnumber+1
    group_frequency(i)=group_frequency(i-1)*exp(dlnnu_group)
end do
write(*,'(5ES15.4E2)') numin,numax
!write(*,'(21ES15.4E2)') group_frequency
call generate_fld_mg(groupnumber,group_frequency)
!do i=1,groupnumber
!    write(*,'(101ES15.4E2)') rad_mg(:,i)
!end do
temp=1d4
do i=1,groupnumber
    call radpower_mg(i,temp,radpower(i))
    call dradpowerdT_mg(i,temp,dradpowerdT(i))
end do
print *,sum(radpower),planck_function(temp)
print *,sum(dradpowerdT),a_rad*c_light/pi*temp**3

contains

subroutine generate_fld_mg(groupnumber,group_frequency)
    !group_frequency gives the lower and upper bound of each rad group
    !each group is further divided into nradcells
    integer :: groupnumber,i,j
    real(8), dimension(:), allocatable :: group_frequency
    real(8) :: dlnnu,multiply
    allocate(rad_mg(nradcell+1,groupnumber))
    do j=1,groupnumber
        dlnnu=log(group_frequency(j+1)/group_frequency(j))/nradcell
        multiply=exp(dlnnu)
        rad_mg(1,j)=group_frequency(j)
        do i=2,nradcell+1
            rad_mg(i,j)=rad_mg(i-1,j)*multiply
        end do
    end do
end subroutine generate_fld_mg

function planck_function_peak_frequency(temp)
    real(8) :: planck_function_peak_frequency,temp
    planck_function_peak_frequency=5.879d10*temp
end function planck_function_peak_frequency

function planck_law_frequency_dlnnu(temp,nu)
    !nu in s^{-1}
    real(8) :: temp,nu,planck_law_frequency_dlnnu,s
    s=h_planck*nu/kb/temp
    if (s<1d-6) then
        planck_law_frequency_dlnnu=2d0*nu**3/c_light**2*kb*temp
    else if (s>20d0) then
        planck_law_frequency_dlnnu=2d0*h_planck*nu**4/c_light**2/exp(s)
    else
        planck_law_frequency_dlnnu=2d0*h_planck*nu**4/c_light**2/(exp(s)-1d0)
    end if
end function planck_law_frequency_dlnnu

function planck_law_dfrequency_dlnnu(temp,nu)
    real(8) :: temp,nu,planck_law_dfrequency_dlnnu,s
    s=h_planck*nu/kb/temp
    if (s<1d-6) then
        planck_law_dfrequency_dlnnu=2d0*nu**3/c_light**2*kb
    else if (s>20d0) then
        planck_law_dfrequency_dlnnu=2d0*h_planck*nu**4/c_light**2/exp(s)*s/temp
    else
        planck_law_dfrequency_dlnnu=2d0*h_planck*nu**4/c_light**2/(exp(s)-1d0)**2*exp(s)*s/temp
    end if
end function planck_law_dfrequency_dlnnu

subroutine erad_mg(igroup,temp,erad)
    integer :: igroup,i
    real(8) :: temp,radpower,nu1,nu2,eradnu1,eradnu2,lnnu1,lnnu2
    erad=0d0
    do i=1,nradcell
        nu1=rad_mg(i,igroup)
        nu2=rad_mg(i+1,igroup)
        lnnu1=log(nu1)
        lnnu2=log(nu2)
        eradnu1=4d0*pi*planck_law_frequency_dlnnu(temp,nu1)/c_light
        eradnu2=4d0*pi*planck_law_frequency_dlnnu(temp,nu2)/c_light
        erad=erad+(eradnu1+eradnu2)/2*(lnnu2-lnnu1)
    end do
end subroutine erad_mg

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

function planck_function(temp)
    !rad intensity integrated over frequency and assume blackbody
    real(8) :: planck_function,temp
    planck_function=a_rad*c_light/4d0/pi*temp**4
end function planck_function

end program multigroup
