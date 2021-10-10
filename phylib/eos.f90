module eos
use datastructure
use mathlib
use phylib
#if ieos==1
    use constant_gamma_eos
#elif ieos==2
    use eos_analytic
#elif ieos==3
    use eos_tabulated
#endif
implicit none

contains


subroutine initialize_eos_environment()
#if     ieos==2
    !eos_hllc_analytic
    call initialize_x_ratio()
    allocate(temp_division(environment%nrho,4))
    call generate_rho_t_eos()
    call generate_temperature_division(environment%rho_eos,environment%t_eos,temp_division)
    if (rank==0) print *,'using analytic eos'
#endif
end subroutine initialize_eos_environment

subroutine convert_u_to_w_block(blk,mark)
    !convert u_field to w_field, all domain
    !given u, calculate w, temp and, egv
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: u(5),w(5),temp,egv
    character(len=128) :: alert
    integer, optional :: mark
    if (nd==1) then
        do i=1,blk_size_nx
#if         ieos==2
            !eos_hllc_analytic
            u=blk%u(i,1,1,1:5)
            call eos_hllc_analytic_utow(u,w,temp,egv)
            blk%w(i,1,1,1:5)=w
            blk%temp(i,1,1)=temp
            blk%egv(i,1,1)=egv
#elif       ieos==1
            u=blk%u(i,1,1,1:5)
            call utow(u,w)
            blk%w(i,1,1,1:5)=w
            blk%temp(i,1,1)=idealgas_temp_p(w(1),w(5))
            blk%egv(i,1,1)=idealgas_egv(w(5))
#endif
        end do
    else if (nd==2) then
        do j=1,blk_size_ny
            do i=1,blk_size_nx
#if             ieos==2
                !eos_hllc_analytic
                u=blk%u(i,j,1,1:5)
                call eos_hllc_analytic_utow(u,w,temp,egv)
                blk%w(i,j,1,1:5)=w
                blk%temp(i,j,1)=temp
                blk%egv(i,j,1)=egv
#elif           ieos==1
                u=blk%u(i,j,1,1:5)
                call utow(u,w)
                blk%w(i,j,1,1:5)=w
                blk%temp(i,j,1)=idealgas_temp_p(w(1),w(5))
                blk%egv(i,j,1)=idealgas_egv(w(5))
#endif
            end do
        end do
    end if
end subroutine convert_u_to_w_block

subroutine convert_w_to_u_block(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: u(5),w(5),temp,egv
    if (nd==1) then
    else if (nd==2) then
        do j=blk_ylb,blk_yub
            do i=blk_xlb,blk_xub
#if             ieos==1
                w=blk%w(i,j,1,1:5)
                call wtou(w,u)
                blk%u(i,j,1,1:5)=u
                blk%temp(i,j,1)=idealgas_temp_p(w(1),w(5))
                blk%egv(i,j,1)=idealgas_egv(w(5))
#elif           ieos==2
#endif
            end do
        end do
    end if
end subroutine convert_w_to_u_block

subroutine convert_u_to_source_cell(u,source)
    real(8), dimension(5) :: u,source
    character(len=128) :: alert
#if     ieos==1
    call utosource(u,source)
#elif   ieos==2
    call eos_analytic_utosource(u,source)
#else
#endif
end subroutine convert_u_to_source_cell

subroutine convert_u_to_source_block(blk)
    !convert blk%u to blk%source
    type(blockdef), pointer :: blk
    real(8) :: u(5),s(5)
    integer :: i,j
    if (nd==1) then
        do i=1,blk_size_nx
            u=blk%u(i,1,1,1:5)
#if     ieos==1
            call utosource(u,s)
#elif   ieos==2
            call eos_analytic_utosource(u,s)
#endif
            blk%source(i,1,1,1:5)=s
        end do
    else if (nd==2) then
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                u=blk%u(i,j,1,1:5)
#if     ieos==1
                call utosource(u,s)
#elif   ieos==2
#endif
                blk%source(i,j,1,1:5)=s
            end do
        end do
    end if
end subroutine convert_u_to_source_block

subroutine convert_source_to_u_block(blk)
    !convert blk%source to blk%u
    type(blockdef), pointer :: blk
    real(8) :: u(5),s(5)
    integer :: i,j
    if (nd==1) then
        do i=1,blk_size_nx
            s=blk%source(i,1,1,1:5)
#if     ieos==1
            call sourcetou(s,u)
#elif   ieos==2
            call eos_analytic_sourcetou(s,u)
#endif
            blk%u(i,1,1,1:5)=u
        end do
    else if (nd==2) then
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                s=blk%source(i,j,1,1:5)
#if     ieos==1
                call sourcetou(s,u)
#elif   ieos==2
#endif
                blk%u(i,j,1,1:5)=u
            end do
        end do
    end if
end subroutine convert_source_to_u_block

subroutine convert_source_to_u_cell(source,u)
    real(8), dimension(5) :: u,source
    character(len=128) :: alert
#if     ieos==1
    call sourcetou(source,u)
#elif   ieos==2
    call eos_analytic_sourcetou(source,u)
#else
#endif
end subroutine convert_source_to_u_cell

subroutine calculate_entropy_block(blk)
    !used only in output step
    !erg/g
    type(blockdef), pointer :: blk
    integer :: i
    real(8) :: rho,temp,ntot,species4(4)
    real(8), allocatable :: species5(:)
    if (nd==1) then
#if     ieos==1
        do i=1,blk_size_nx
            rho=blk%w(i,1,1,1)
            temp=blk%temp(i,1,1)
            blk%entropy(i,1,1)=idealgas_entropy(rho,temp)
        end do
#elif   ieos==2
        do i=1,blk_size_nx
            rho=blk%w(i,1,1,1)
            temp=blk%temp(i,1,1)
            blk%entropy(i,1,1)=srhot(rho,temp)
            !blk%entropy(i,1,1)=srhot_per_particle(rho,temp)
        end do
#endif
    else if (nd==2) then
    end if
end subroutine calculate_entropy_block

subroutine get_temp_from_egv_rho(egv,w,temp)
    real(8), allocatable, dimension(:,:,:) :: egv,temp
    real(8), allocatable, dimension(:,:,:,:) :: w
end subroutine get_temp_from_egv_rho

subroutine convert_rho_temp_to_u_w_block(blk)
    type(blockdef), pointer :: blk
    real(8), allocatable, dimension(:,:,:) :: egv_old
    integer :: i,j
#if     ieos==1
    call allocate_cell_data_block(egv_old)
    egv_old=blk%egv
    if (nd==1) then
        do i=1,blk_size_nx
            blk%egv(i,1,1)=blk%w(i,1,1,1)*blk%temp(i,1,1)*kb/(gamma_gas-1)/maw/amu
            blk%u(i,1,1,5)=blk%u(i,1,1,5)-egv_old(i,1,1)+blk%egv(i,1,1)
            blk%w(i,1,1,5)=calculate_p_from_egv_rho(blk%egv(i,1,1),blk%w(i,1,1,1))
        end do
    else if (nd==2) then
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                blk%egv(i,j,1)=blk%w(i,j,1,1)*blk%temp(i,j,1)*kb/(gamma_gas-1)/maw/amu
                blk%u(i,j,1,5)=blk%u(i,j,1,5)-egv_old(i,j,1)+blk%egv(i,j,1)
                blk%w(i,j,1,5)=calculate_p_from_egv_rho(blk%egv(i,j,1),blk%w(i,j,1,1))
            end do
        end do
    end if
    deallocate(egv_old)
#elif   ieos==2
    call allocate_cell_data_block(egv_old)
    egv_old=blk%egv
    if (nd==1) then
        do i=1,blk_size_nx
            blk%egv(i,1,1)=egvrhot(blk%w(i,1,1,1),blk%temp(i,1,1))
            blk%u(i,1,1,5)=blk%u(i,1,1,5)-egv_old(i,1,1)+blk%egv(i,1,1)
            blk%w(i,1,1,5)=calculate_p_from_egv_rho(blk%egv(i,1,1),blk%w(i,1,1,1))
        end do
    else if (nd==2) then
    end if
    deallocate(egv_old)
#endif
end subroutine convert_rho_temp_to_u_w_block

subroutine get_cv_block(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
#if     ieos==1
    if (nd==1) then
        do i=1,blk_size_nx
            blk%cv(i,1,1)=blk%w(i,1,1,1)*kb/(gamma_gas-1)/maw/amu
        end do
    else if (nd==2) then
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                blk%cv(i,j,1)=blk%w(i,j,1,1)*kb/(gamma_gas-1)/maw/amu
            end do
        end do
    end if
#elif   ieos==2
    if (nd==1) then
        do i=1,blk_size_nx
            blk%cv(i,1,1)=eos_h_cv(blk%w(i,1,1,1),blk%temp(i,1,1))
        end do
    else if (nd==2) then
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                blk%cv(i,j,1)=eos_h_cv(blk%w(i,1,1,1),blk%temp(i,1,1))
            end do
        end do
    end if
#endif
end subroutine get_cv_block

function calculate_temp_from_p_rho(p,rho)
    real(8) :: calculate_temp_from_p_rho,p,rho
#if     ieos==1
    calculate_temp_from_p_rho=idealgas_temp_p(rho,p)
#elif   ieos==2
    calculate_temp_from_p_rho=solvetp(p,rho)
#endif
end function calculate_temp_from_p_rho

function calculate_temp_from_egv_rho(egv,rho)
    real(8) :: egv,rho,calculate_temp_from_egv_rho
#if     ieos==1
    calculate_temp_from_egv_rho=idealgas_temp_egv(rho,egv)
#elif   ieos==2
    calculate_temp_from_egv_rho=solvetegv(egv,rho)
#endif
end function calculate_temp_from_egv_rho

function calculate_egv_from_rho_temp(rho,temp)
    real(8) :: calculate_egv_from_rho_temp,rho,temp
#if     ieos==1
    calculate_egv_from_rho_temp=rho*temp*kb/(gamma_gas-1)/maw/amu
#elif   ieos==2
#endif
end function calculate_egv_from_rho_temp

function calculate_p_from_egv_rho(egv,rho)
    real(8) :: egv,rho,calculate_p_from_egv_rho,temp
#if     ieos==1
    calculate_p_from_egv_rho=idealgas_p(egv)
#elif   ieos==2
    temp=solvetegv(egv,rho)
    calculate_p_from_egv_rho=prhot(rho,temp)
#endif
end function calculate_p_from_egv_rho

function calculate_p_from_rho_t(rho,t)
    real(8) :: rho,t,calculate_p_from_rho_t
#if     ieos==1
    calculate_p_from_rho_t=idealgas_prhot(rho,t) 
#elif   ieos==2
    calculate_p_from_rho_t=prhot(rho,t)
#endif
end function calculate_p_from_rho_t

function calculate_adiabatic_cs(rho,t)
    real(8) :: rho,t,calculate_adiabatic_cs
#if     ieos==1
    calculate_adiabatic_cs=sqrt(gamma_gas*kb*t/maw/amu)
#elif   ieos==2
    calculate_adiabatic_cs=adiabatic_cs(rho,t)
#elif   ieos==3
#else
    calculate_adiabatic_cs=0d0
#endif
end function calculate_adiabatic_cs

subroutine convert_u_to_w_cell(u,w,temp,egv)
    real(8) :: w(5),u(5),temp,egv
    character(len=128) :: alert
#if     ieos==1
    call utow(u,w)
    temp=idealgas_temp_p(w(1),w(5))
    egv=idealgas_egv(w(5))
#elif   ieos==2
    call eos_hllc_analytic_utow(u,w,temp,egv)
#elif   ieos==3
    call eos_utow(u,w)
#else
    alert='ieos number wrong'
    call abort_achilles(alert)
#endif
end subroutine convert_u_to_w_cell

subroutine examine_temperature(temp)
    real(8), dimension(:,:,:), allocatable :: temp
    integer :: i,j,k
    if (nd==1) then
        do i=xlb,xub
            if (temp(i,1,1)<=0d0) then
                print *,'negative temperature found at index', i
                stop
            end if
        end do
    else if (nd==2) then
    else
    end if
end subroutine examine_temperature


#if     ieos==3
!subroutine cal_species_h_he(info)
!    !calculate all species (H2,HI,HII,HeI,HeII,HeIII,e-)
!    type(infodef) :: info
!    real(8) :: x,y,mu,species(7),mu1,mu2,mu3,mu4,mu5
!    integer :: i,j,k,chunk
!    x=info%h_ratio
!    y=info%he_ratio
!    mu1=1d0/(x/muh2+y/muhe)
!    mu2=1d0/(x/muh+y/muhe)
!    mu3=1d0/(2d0*x/muh+y/muhe)
!    mu4=1d0/(2d0*x/muh+2d0*y/muhe)
!    mu5=1d0/(2d0*x/muh+3d0*y/muhe)
!    do i=1,nx
!        do j=1,ny
!            do k=1,nz
!                info%mu(i,j,k)=eos_tabulated_maw_2d(info%w(i,j,k,1),info%temp(i,j,k))
!                mu=info%mu(i,j,k)
!                if (mu.le.mu1*1.0001.and.mu.gt.mu1) then
!                    call h2_hI_heI(x,y,mu1,species)
!                else if (mu.le.mu1.and.mu.gt.mu2) then
!                    call h2_hI_heI(x,y,mu,species)
!                else if (mu.le.mu2.and.mu.gt.mu3) then
!                    call hI_hII_heI_e(x,y,mu,species)
!                else if (mu.le.mu3.and.mu.gt.mu4) then
!                    call hII_heI_heII_e(x,y,mu,species)
!                else if (mu.le.mu4.and.mu.gt.mu5) then
!                    call hII_heII_heIII_e(x,y,mu,species)
!                else if (mu.le.mu5.and.mu.gt.mu5/1.0001) then
!                    call hII_heII_heIII_e(x,y,mu5,species)
!                else
!                    print *, 'mu out of bound',mu,mu1,mu2,mu3,mu4,mu5,x,y
!                    stop
!                end if
!                info%species(i,j,k,:)=species
!            end do
!        end do
!    end do
!end subroutine cal_species_h_he
#endif

end module eos
