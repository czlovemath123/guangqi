module phylib
use datastructure
use mathlib
implicit none

!**************interp method, EOS Riemann solver related tables and quantities***********
!type(table2d), public :: cs_trho,eg_trho,t_prho,p_trho,gamma_trho,eg_rhop,rho_egp,maw_trho,t_egrho
!type(table1d), public :: gamma_t,egm_t,g_t,t_g,cs_t,maw_t,g_egm,egm_g,t_egm,cs_g
!real(8), public :: rho_min_bound,rho_max_bound,p_min_bound,p_max_bound,cs_min_bound,cs_max_bound,  &
!    t_min_bound,t_max_bound,gamma_max,gamma_min,g_min_bound,g_max_bound
!**************interp method, EOS Riemann solver related tables and quantities***********

real(8) :: dust_beta,opacity_gas,opacity_dust,dens_threshold_opacity,dens_cooling_threshold,  &
    dust_formation_length_scale,I_dust_threshold,r_dust_formation,dust_radcool(4),dust_gas_ratio
integer :: idustform
public :: dust_beta,opacity_gas,opacity_dust,dens_threshold_opacity,dens_cooling_threshold,  &
    dust_formation_length_scale,I_dust_threshold,r_dust_formation,dust_radcool,idustform,dust_gas_ratio

!**************radiation transfer related tables and quantities**************************
type(table1d), public :: NK_h2o_l0
type(table2d), public :: NK_h2o_llte,NK_h2o_n0_5,NK_h2o_alpha
type(table2d), public :: kappa_lowT,kappa_highT
type(table2d), public :: rosseland_gas_opacity_table,planck_gas_opacity_table
type(table2d), public :: rosseland_dust_opacity_table,planck_dust_opacity_table
real(8), public :: critical_dtao
real(8), public :: t_min_kappa,t_max_kappa,r_min_kappa,r_max_kappa  !t=log10, r=log10(rho)-3log10(t)+18
!**************radiation transfer related tables and quantities**************************

contains

!*************boundary specification. the quantities needed for specification include****
!*********w(1:5), temp, egv, and u(1:5) at the boundary for all ieos*********************
!*********object oriented subroutines. all boundaries are objects************************
    
subroutine apply_hydro_condition(sub)
    !assign value to all cells including guard cells in all blocks
    type(blockdef), pointer :: blk
    procedure(condition_hydro) :: sub
    real(8) :: pos(3),t,w(5),u(5),temp,egv,Frad(3),feddington
    real(8), dimension(:), allocatable :: Erad
    integer :: ijk(3),i,j,k,iblk
    blk=>blk_processor_head
    t=time_sys%t
    if (nd==1) then
        do iblk=1,np_nblk(rank+1)
            do i=blk_xlb,blk_xub
                pos=(/blk%x_center(i),0d0,0d0/)
                call sub(pos,t,w,u,temp,egv)
                blk%w(i,1,1,1:5)=w
                blk%u(i,1,1,1:5)=u
                blk%temp(i,1,1)=temp
                blk%egv(i,1,1)=egv
            end do
            blk=>blk%blk_next
        end do
    else if (nd==2) then
        do iblk=1,np_nblk(rank+1)
            do j=blk_ylb,blk_yub
                do i=blk_xlb,blk_xub
                    pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                    call sub(pos,t,w,u,temp,egv)
                    blk%w(i,j,1,1:5)=w
                    blk%u(i,j,1,1:5)=u
                    blk%temp(i,j,1)=temp
                    blk%egv(i,j,1)=egv
                end do
            end do
            if (iblk/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
    else
    end if
    nullify(blk)
end subroutine apply_hydro_condition

subroutine apply_rad_condition(sub)
    type(blockdef), pointer :: blk
    procedure(condition_rad) :: sub
    real(8) :: pos(3),t,Frad(3),feddington,temp
    real(8), dimension(:), allocatable :: Erad
    integer :: ijk(3),i,j,k,iblk
    blk=>blk_processor_head
    t=time_sys%t
    if (nd==1) then
        do iblk=1,np_nblk(rank+1)
            do i=blk_xlb,blk_xub
                pos=(/blk%x_center(i),0d0,0d0/)
                if (iradiation==1) then
                    allocate(Erad(1))
                    call sub(pos,t,Erad,Frad,feddington)
                    blk%Erad(i,1,1)=Erad(1)
                    deallocate(Erad)
                else if (iradiation==3) then
                else if (iradiation==4) then
                    allocate(Erad(1))
                    call sub(pos,t,Erad,Frad,feddington)
                    blk%Erad(i,1,1)=Erad(1)
                    deallocate(Erad)
                else if (iradiation==5) then
                end if
            end do
            blk=>blk%blk_next
        end do
    else if (nd==2) then
        do iblk=1,np_nblk(rank+1)
            do j=blk_ylb,blk_yub
                do i=blk_xlb,blk_xub
                    pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                    if (iradiation==1) then
                        allocate(Erad(1))
                        call sub(pos,t,Erad,Frad,feddington)
                        blk%Erad(i,j,1)=Erad(1)
                        deallocate(Erad)
                    else if (iradiation==3) then
                    else if (iradiation==4) then
                        allocate(Erad(1))
                        temp=blk%temp(i,j,1)
                        call sub(pos,t,Erad,Frad,feddington)
                        blk%Erad(i,j,1)=a_rad*temp**4d0
                        deallocate(Erad)
                    end if
                end do
            end do
            if (iblk/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
    else
    end if
    nullify(blk)
end subroutine apply_rad_condition

subroutine apply_scalar_condition(sub,sub1)
    !assign scalar value to all domain
    type(blockdef), pointer :: blk
    procedure(condition_scalar), pointer :: sub
    procedure(extractarray) :: sub1
    real(8) :: pos(3),t,scalar
    real(8), dimension(:,:,:), pointer :: ap
    integer :: ijk(3),i,j,k,iblk
    blk=>blk_processor_head
    if (nd==1) then
    else if (nd==2) then
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        do iblk=1,np_nblk(rank+1)
            call sub1(blk,ap)
            do j=blk_ylb,blk_yub
                do i=blk_xlb,blk_xub
                    pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                    call sub(pos,t,scalar)
                    ap(i,j,1)=scalar
                end do
            end do
            if (iblk/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
        nullify(ap)
    end if
end subroutine apply_scalar_condition

subroutine apply_condition_to_domain(sub,subrad)
    type(blockdef), pointer :: blk
    procedure(condition_domain) :: sub
    procedure(condition_domainrad), optional :: subrad
    integer :: i,j,k,ijk(3),iblk
    real(8) :: pos(3),w(5),u(5),temp,egv,t,Erad,Frad(3),feddington
    logical :: assign_value
    character(len=128) :: alert
    blk=>blk_processor_head
    t=time_sys%t
    if (nd==1) then
        do iblk=1,np_nblk(rank+1)
            do i=blk_xlb,blk_xub
                pos=(/blk%x_center(i),0d0,0d0/)
                call sub(pos,t,w,u,temp,egv,assign_value)
                if (assign_value) then
                    blk%w(i,1,1,1:5)=w
                    blk%u(i,1,1,1:5)=u
                    blk%temp(i,1,1)=temp
                    blk%egv(i,1,1)=egv
                end if
                if (iradiation==3) then
                    call subrad(pos,t,Erad,Frad,feddington,assign_value)
                    if (assign_value) then
                        blk%Erad(i,1,1)=Erad
                        blk%Fradx(i,1,1)=Frad(1)
                        blk%feddington(i,1,1)=feddington
                    end if
                end if
            end do
            if (iblk/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
    else if (nd==2) then
        do iblk=1,np_nblk(rank+1)
            do j=blk_ylb,blk_yub
                do i=blk_xlb,blk_xub
                    pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                    call sub(pos,t,w,u,temp,egv,assign_value)
                    if (assign_value) then
                        blk%w(i,j,1,1:5)=w
                        blk%u(i,j,1,1:5)=u
                        blk%temp(i,j,1)=temp
                        blk%egv(i,j,1)=egv
                    end if
                    if (iradiation==3) then
                        call subrad(pos,t,Erad,Frad,feddington,assign_value)
                        if (assign_value) then
                            blk%Erad(i,j,1)=Erad
                            blk%Fradx(i,j,1)=Frad(1)
                            blk%Frady(i,j,1)=Frad(2)
                            blk%feddington(i,j,1)=feddington
                        end if
                    end if
                end do
            end do
            if (iblk/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
    else
    end if
end subroutine apply_condition_to_domain

subroutine cal_species_h_he_indivi(x,y,mu,species)
    real(8) :: x,y,mu,mu1,mu2,mu3,mu4,mu5,species(7)
    mu1=1d0/(x/muh2+y/muhe)
    mu2=1d0/(x/muh+y/muhe)
    mu3=1d0/(2d0*x/muh+y/muhe)
    mu4=1d0/(2d0*x/muh+2d0*y/muhe)
    mu5=1d0/(2d0*x/muh+3d0*y/muhe)
    if (mu.le.mu1*1.0001.and.mu.gt.mu1) then
        call h2_hI_heI(x,y,mu1,species)
    else if (mu.le.mu1.and.mu.gt.mu2) then
        call h2_hI_heI(x,y,mu,species)
    else if (mu.le.mu2.and.mu.gt.mu3) then
        call hI_hII_heI_e(x,y,mu,species)
    else if (mu.le.mu3.and.mu.gt.mu4) then
        call hII_heI_heII_e(x,y,mu,species)
    else if (mu.le.mu4.and.mu.gt.mu5) then
        call hII_heII_heIII_e(x,y,mu,species)
    else if (mu.le.mu5.and.mu.gt.mu5/1.0001) then
        call hII_heII_heIII_e(x,y,mu5,species)
    else
        print *, 'mu out of bound',mu,mu1,mu2,mu3,mu4,mu5,x,y
        stop
    end if
end subroutine cal_species_h_he_indivi

subroutine h2_hI_heI(x,y,mu,species)
    !h2 disassociates to become hI
    real(8) :: x,y,a,mu,species(7)
    a=(1-mu*x/muh2-mu*y/muhe)/(mu*x*(1d0/muh-1d0/muh2))
    species(1)=(1d0-a)*x
    species(2)=a*x
    species(3)=0d0
    species(4)=y
    species(5)=0d0
    species(6)=0d0
    species(7)=0d0
end subroutine h2_hI_heI

subroutine hI_hII_heI_e(x,y,mu,species)
    !hI becomes hII and e-
    real(8) :: x,y,a,mu,species(7)
    a=(1-mu*x/muh-mu*y/muhe)/(mu*x*(2d0/muh-1d0/muh))
    species(1)=0d0
    species(2)=(1d0-a)*x
    species(3)=a*x
    species(4)=y
    species(5)=0d0
    species(6)=0d0
    species(7)=a*x
end subroutine hI_hII_heI_e

subroutine hII_heI_heII_e(x,y,mu,species)
    !heI becomes heII and e-
    real(8) :: x,y,a,mu,species(7)
    a=(1-2*mu*x/muh-mu*y/muhe)/(mu*y*(2d0/muhe-1d0/muhe))
    species(1)=0d0
    species(2)=0d0
    species(3)=x
    species(4)=(1d0-a)*y
    species(5)=a*y
    species(6)=0d0
    species(7)=x+a*y
end subroutine hII_heI_heII_e

subroutine hII_heII_heIII_e(x,y,mu,species)
    !heII becomes heIII and e-
    real(8) :: x,y,a,mu,species(7)
    a=(1d0-2d0*mu*x/muh-2d0*mu*y/muhe)/(mu*y*(3d0/muhe-2d0/muhe))
    species(1)=0d0
    species(2)=0d0
    species(3)=x
    species(4)=0d0
    species(5)=(1d0-a)*y
    species(6)=a*y
    species(7)=x+2d0*a*y+(1d0-a)*y
end subroutine hII_heII_heIII_e


!******************************************************************************
!*******************************end EOS part***********************************
!******************************************************************************

subroutine initialize_star(star,m,temp,r,period,spin)
    type(stellar) :: star
    real(8) :: m,loc(3),temp,r
    real(8), optional :: period,spin
    loc=0d0
    star%core%mass=m
    star%core%xyz=loc
    star%teff=temp
    star%radius=r
    star%luminosity=4*pi*r*r*sigma_sb*temp**4
    star%lum_atm=0d0
    if (present(period)) then
        star%period=period   !if the star has some periodic property
    else
        star%period=t_hubble
    end if
    if (present(spin)) then
        star%spin=spin
    else
        star%spin=0d0
    end if
end subroutine initialize_star

function blackbody_rad_power(temp)
    !Stefan-Boltzmann law, the intensity has a 1/pi relation
    real(8) :: blackbody_rad_power,temp
    blackbody_rad_power=sigma_sb*temp**4
end function blackbody_rad_power

function planck_function(temp)
    !rad intensity integrated over frequency and assume blackbody
    real(8) :: planck_function,temp
    planck_function=a_rad*c_light/4d0/pi*temp**4
end function planck_function

function planck_law_wavelength(temp,x)
    !x in cm
    real(8) :: temp,x,planck_law_wavelength
    planck_law_wavelength=2*h_planck*c_light**2/x**5/(exp(h_planck*c_light/x/kb/temp)-1)
end function planck_law_wavelength

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

function rad_energy_density(temp)
    !radiation energy density assume blackbody
    real(8) :: rad_energy_density,temp
    rad_energy_density=a_rad*temp**4
end function rad_energy_density

subroutine rad_energy_density_mg_all(Erad_mg,temp)
    real(8) :: temp
    real(8), dimension(:), allocatable :: Erad_mg
    integer :: i
    do i=1,nmg
        Erad_mg(i)=rad_energy_density_mg(i,temp)
    end do
end subroutine rad_energy_density_mg_all

function rad_energy_density_mg(igroup,temp)
    integer :: igroup,i
    real(8) :: temp,rad_energy_density_mg,nu1,nu2,eradnu1,eradnu2,lnnu1,lnnu2
    rad_energy_density_mg=0d0
    do i=1,nradcell
        nu1=rad_mg(i,igroup)
        nu2=rad_mg(i+1,igroup)
        lnnu1=log(nu1)
        lnnu2=log(nu2)
        eradnu1=4d0*pi*planck_law_frequency_dlnnu(temp,nu1)/c_light
        eradnu2=4d0*pi*planck_law_frequency_dlnnu(temp,nu2)/c_light
        rad_energy_density_mg=rad_energy_density_mg+(eradnu1+eradnu2)/2*(lnnu2-lnnu1)
    end do
end function rad_energy_density_mg

subroutine block_environment_setup()
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: m,r,pos(3),dxyz(3)
end subroutine block_environment_setup

function sum_eg()
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: sum_eg
    sum_eg=0d0
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            sum_eg=sum_eg+blk%vol(j,1,1)*blk%egv(j,1,1)
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end function sum_eg

function sum_ek()
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: sum_ek
    sum_ek=0d0
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            sum_ek=sum_ek+blk%vol(j,1,1)*(blk%u(j,1,1,5)-blk%egv(j,1,1))
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end function sum_ek

function sum_er()
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: sum_er
    sum_er=0d0
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            sum_er=sum_er+blk%vol(j,1,1)*blk%Erad(j,1,1)
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end function sum_er

function sum_gpotential()
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: sum_gpotential
    sum_gpotential=0d0
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        do j=1,blk_size_nx
            sum_gpotential=sum_gpotential+blk%vol(j,1,1)*blk%w(j,1,1,1)*blk%gpotential(j,1,1)
        end do
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end function sum_gpotential

subroutine spherical2d_vphi(blk,vphi)
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: vphi
    vphi=>blk%vphi
end subroutine spherical2d_vphi

end module phylib
