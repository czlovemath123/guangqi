module eos_hllc_roe_analytic
use datastructure
use phylib
use mathlib
use eos_analytic
implicit none
contains

!P. Batten, N. Clarke, C. Lambert, and D. M. Causon "On the choince of wavespeeds for the HLLC Riemann solver" 1997

subroutine eos_hllc_roe_analytic_roe_average(wl,wr,egv,q_roe,c_roe)
    !use Roe average Riemann solver to find q_roe, c_roe
    real(8), dimension(5) :: wl,wr,w_roe
    real(8) :: egv(2),hl,hr,rhot,ut,vt,wt,ht,h,tt,q_roe,c_roe
    real(8) :: sqrtrhol,sqrtrhor,sqrtrho,hv
    hl=(half*wl(1)*(wl(2)**2+wl(3)**2+wl(4)**2)+egv(1)+wl(5))/wl(1)
    hr=(half*wr(1)*(wr(2)**2+wr(3)**2+wr(4)**2)+egv(2)+wr(5))/wr(1)
    rhot=sqrt(wl(1)*wr(1))
    sqrtrhol=sqrt(wl(1))
    sqrtrhor=sqrt(wr(1))
    sqrtrho=sqrt(wl(1))+sqrt(wr(1))
    ut=(sqrtrhol*wl(2)+sqrtrhor*wr(2))/sqrtrho
    vt=(sqrtrhol*wl(3)+sqrtrhor*wr(3))/sqrtrho
    wt=(sqrtrhol*wl(4)+sqrtrhor*wr(4))/sqrtrho
    ht=(sqrtrhol*hl+sqrtrhor*hr)/sqrtrho
    hv=(ht-half*(ut**2+vt**2+wt**2))*rhot
    tt=solveth(hv,rhot)
    c_roe=adiabatic_cs(rhot,tt)
    q_roe=ut
end subroutine eos_hllc_roe_analytic_roe_average

subroutine eos_hllc_roe_analytic_flux(wl,wr,temp,egv,flux,vlocalmax)
    !extended infomation:rho,vx,vy,vz,p,temp,egv
    !given left and right states (wl,wr), temp, and egv, calculate flux and vlocalmax
    real(8), dimension(5) :: wl,wr,flux
    real(8) :: sl,sr,sm,vlocalmax,temp(2),egv(2)
    real(8) :: ql,qr,q_roe,cl,cr,c_roe,dsl,dsr,dsm
    real(8) :: rhom,pm,momum,momvm,momwm,el,er,em
    call eos_hllc_roe_analytic_roe_average(wl,wr,egv,q_roe,c_roe)
    ql=wl(2)
    qr=wr(2)
    cl=adiabatic_cs(wl(1),temp(1))
    cr=adiabatic_cs(wr(1),temp(2))
    sl=min(ql-cl,q_roe-c_roe)
    sr=max(qr+cr,q_roe+c_roe)
    dsl=sl-ql
    dsr=sr-qr
    sm=(wr(1)*qr*dsr-wl(1)*ql*dsl+wl(5)-wr(5))/(wr(1)*dsr-wl(1)*dsl)
    if (sl>zero) then
        call eos_hllc_analytic_wtoxflux(wl,egv(1),flux)
    else if (sl<=zero.and.sm>zero) then
        dsm=sl-sm
        rhom=wl(1)*dsl/dsm
        pm=wl(5)-wl(1)*dsl*(ql-sm)
        momum=(dsl*wl(1)*wl(2)+(pm-wl(5)))/dsm
        momvm=dsl*wl(1)*wl(3)/dsm
        momwm=dsl*wl(1)*wl(4)/dsm
        el=egv(1)+half*wl(1)*(wl(2)**2+wl(3)**2+wl(4)**2)
        em=(dsl*el-ql*wl(5)+pm*sm)/dsm
        flux(1)=rhom*sm
        flux(2)=momum*sm+pm
        flux(3)=momvm*sm
        flux(4)=momwm*sm
        flux(5)=(em+pm)*sm
    else if (sm<=zero.and.sr>zero) then
        dsm=sr-sm
        rhom=wr(1)*dsr/dsm
        pm=wr(5)-wr(1)*dsr*(qr-sm)
        momum=(dsr*wr(1)*wr(2)+(pm-wr(5)))/dsm
        momvm=dsr*wr(1)*wr(3)/dsm
        momwm=dsr*wr(1)*wr(4)/dsm
        er=egv(2)+half*wr(1)*(wr(2)**2+wr(3)**2+wr(4)**2)
        em=(dsr*er-qr*wr(5)+pm*sm)/dsm
        flux(1)=rhom*sm
        flux(2)=momum*sm+pm
        flux(3)=momvm*sm
        flux(4)=momwm*sm
        flux(5)=(em+pm)*sm
    else
        call eos_hllc_analytic_wtoxflux(wr,egv(2),flux)
    end if
    vlocalmax=max(abs(sl),abs(sr))
end subroutine eos_hllc_roe_analytic_flux

subroutine eos_hllc_roe_analytic_sweep(w,temp,egv,dudt,vglobalmax,dx,dir)
    !given w, temp, egv, and dir, calculate dudt and vglobalmax
    !u is the conserved quantities in dudt
    real(8), dimension(:,:,:,:), allocatable :: w
    real(8), dimension(:,:,:,:), allocatable :: dudt
    real(8), allocatable :: vmaxarray(:,:,:),fluxarray(:,:,:,:),temp(:,:,:),egv(:,:,:)
    real(8), dimension(5) :: wl,wr,flux
    real(8) :: vlocalmax,vglobalmax(nd),dx,temp_local(2),egv_local(2)
    character(len=1), optional :: dir
    character(len=128) :: alert
    integer :: i,j
    if (nd==1) then
        !for 1d, the default is in x direction
        allocate(vmaxarray(nx+1,ny,1),fluxarray(nx+1,ny,1,5))
        vmaxarray=0d0
        fluxarray=0d0
        dudt=0d0
        do i=0,nx
            wl=w(i,1,1,:)
            wr=w(i+1,1,1,:)
            temp_local=(/temp(i,1,1),temp(i+1,1,1)/)
            egv_local=(/egv(i,1,1),egv(i+1,1,1)/)
            !print *,i,temp_local,wl(1),wl(5),wr(1),wr(5)
            call eos_hllc_roe_analytic_flux(wl,wr,temp_local,egv_local,flux,vlocalmax)
            fluxarray(i+1,1,1,:)=flux
            vmaxarray(i+1,1,1)=vlocalmax
        end do
        vglobalmax=maxval(vmaxarray)
        do i=1,nx
            dudt(i,1,1,1:5)=(fluxarray(i,1,1,1:5)-fluxarray(i+1,1,1,1:5))/dx
        end do
        deallocate(vmaxarray,fluxarray)
    else if (nd==2) then
        !for 2d, the default is in x and y direction
        if (dir=='x') then
            allocate(vmaxarray(nx+1,0:ny+1,1),fluxarray(nx+1,0:ny+1,1,5))
            vmaxarray=0d0
            fluxarray=0d0
            dudt=0d0
            do i=0,nx
                do j=0,ny+1
                    wl=w(i,j,1,:)
                    wr=w(i+1,j,1,:)
                    temp_local=(/temp(i,j,1),temp(i+1,j,1)/)
                    egv_local=(/egv(i,j,1),egv(i+1,j,1)/)
                    call eos_hllc_roe_analytic_flux(wl,wr,temp_local,egv_local,flux,vlocalmax)
                    fluxarray(i+1,j,1,:)=flux
                    vmaxarray(i+1,j,1)=vlocalmax
                end do
            end do
            do i=1,nx
                do j=0,ny+1
                    dudt(i,j,1,1:5)=(fluxarray(i,j,1,1:5)-fluxarray(i+1,j,1,1:5))/dx
                end do
            end do
            vglobalmax(1)=maxval(vmaxarray)
            deallocate(vmaxarray,fluxarray)
        else if (dir=='y') then
            allocate(vmaxarray(0:nx+1,ny+1,1),fluxarray(0:nx+1,ny+1,1,5))
            vmaxarray=0d0
            fluxarray=0d0
            dudt=0d0
            do i=0,nx+1
                do j=0,ny
                    wl=w(i,j,1,:)
                    wr=w(i,j+1,1,:)
                    temp_local=(/temp(i,j,1),temp(i,j+1,1)/)
                    egv_local=(/egv(i,j,1),egv(i,j+1,1)/)
                    call rotate_xyz_to_yzx(wl)
                    call rotate_xyz_to_yzx(wr)
                    call eos_hllc_roe_analytic_flux(wl,wr,temp_local,egv_local,flux,vlocalmax)
                    call rotate_yzx_to_xyz(flux)
                    fluxarray(i,j+1,1,:)=flux
                    vmaxarray(i,j+1,1)=vlocalmax
                end do
            end do
            vglobalmax(2)=maxval(vmaxarray)
            do i=0,nx
                do j=1,ny
                    dudt(i,j,1,1:5)=(fluxarray(i,j,1,1:5)-fluxarray(i,j+1,1,1:5))/dx
                end do
            end do
            deallocate(vmaxarray,fluxarray)
        else
            alert='hydro_sweep no such direction'
            call abort_achilles(alert)
        end if
    else
    end if
end subroutine eos_hllc_roe_analytic_sweep

end module eos_hllc_roe_analytic
