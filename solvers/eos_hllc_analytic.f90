module eos_hllc_analytic
use datastructure
use phylib
use mathlib
use eos_analytic
implicit none
contains

function eos_analytic_gleft(wl,p,gamm)
    real(8), dimension(5) :: wl
    real(8) :: p,al,bl,eos_analytic_gleft,gamm
    al=2/(gamm+1)/wl(1)
    bl=(gamm-1)/(gamm+1)*wl(5)
    eos_analytic_gleft=sqrt(al/(p+bl))
end function eos_analytic_gleft

function eos_analytic_gright(wr,p,gamm)
    real(8), dimension(5) :: wr
    real(8) :: p,ar,br,eos_analytic_gright,gamm
    ar=2/(gamm+1)/wr(1)
    br=(gamm-1)/(gamm+1)*wr(5)
    eos_analytic_gright=sqrt(ar/(p+br))
end function eos_analytic_gright

subroutine eos_analytic_hllc_pvrs(wl,wr,pvrs,gammal,gammar)
    real(8), dimension(5) :: wl,wr
    real(8) :: ppvrs,pvrs,rhoa,aa,al,ar,gammal,gammar
    al=sqrt(gammal*wl(5)/wl(1))
    ar=sqrt(gammar*wr(5)/wr(1))
    aa=half*(al+ar)
    rhoa=half*(WL(1)+WR(1))
    ppvrs=half*(WL(5)+WR(5))-half*(WR(2)-WL(2))*rhoa*aa
    pvrs=max(zero,ppvrs)
end subroutine eos_analytic_hllc_pvrs

subroutine eos_analytic_hllc_trrs(wl,wr,p,gammal,gammar,gamm)
    real(8), dimension(5) :: wl,wr
    real(8) :: p,z,al,ar,gammal,gammar,gamm,zeta
    al=sqrt(gammal*wl(5)/wl(1))
    ar=sqrt(gammar*wr(5)/wr(1))
    z=(gamm-1)/two/gamm
    zeta=max((al+ar-gamm*z*(wr(2)-wl(2)))/(al/pow(wl(5),z)+ar/pow(wr(5),z)),0d0)
    p=max(pow(zeta,1d0/z),pfloor)
end subroutine eos_analytic_hllc_trrs

subroutine eos_analytic_hllc_tsrs(wl,wr,p,gammal,gammar,gamm)
    real(8), dimension(5) :: wl,wr
    real(8) :: p,gl,gr,p0,pvrs,gammal,gammar,gamm,x
    call eos_analytic_hllc_pvrs(wl,wr,pvrs,gammal,gammar)
    p0=max(zero,pvrs)
    !x=max(gammal,gammar)
    !gl=eos_analytic_gleft(wl,p0,x)
    !gr=eos_analytic_gright(wr,p0,x)
    gl=eos_analytic_gleft(wl,p0,gammal)
    gr=eos_analytic_gright(wr,p0,gammar)
    p=(gl*wl(5)+gr*wr(5)-(wr(2)-wl(2)))/(gl+gr)
end subroutine eos_analytic_hllc_tsrs

!eos_analytic_hllc_pstar should do the same thing as hllc_pstar
subroutine eos_analytic_hllc_pstar(wl,wr,pstar,gammal,gammar,gamm)
    real(8), dimension(5) :: wl,wr
    real(8) :: pstar,q,qmax,pvrs,trrs,tsrs,pmax,pmin,gammal,gammar,gamm
    call eos_analytic_hllc_pvrs(wl,wr,pvrs,gammal,gammar)
    pmax=max(wl(5),wr(5))
    pmin=min(wl(5),wr(5))
    q=pmax/pmin
    qmax=2.0
    if (q<=qmax.and.pvrs<=pmax.and.pvrs>=pmin) then
        pstar=pvrs
    else
        if (pvrs<pmin) then
            call eos_analytic_hllc_trrs(wl,wr,pstar,gammal,gammar,gamm)
        else
            call eos_analytic_hllc_tsrs(wl,wr,pstar,gammal,gammar,gamm)
        end if
    end if
end subroutine eos_analytic_hllc_pstar

subroutine eos_hllc_analytic_slsrsstar(wl,wr,sl,sr,sstar,egv,gamm,gammal,gammar,pstar)
    !Hybrid estimate in "Restoration of the contact surface in the HLL-Riemann solver"
    !calculate sl, sr, and sstar
    real(8), dimension(5) :: wl,wr
    real(8) :: sl,sr,sl_temp,sr_temp,sstar,pstar,wsl,wsr,al,ar,ql,qr,gamm,gammal,gammar
    real(8) :: sstar_real(2),p_jump,temp(2),egv(2)
    call eos_analytic_hllc_pstar(wl,wr,pstar,gammal,gammar,gamm)
    !call eos_analytic_roe_averaged_rs(wl,wr,egv,pstar,gamm)
    !calculate sl,sr and sstar
    if (pstar>wl(5)) then
        !left shock
        al=sqrt(gammal*wl(5)/wl(1))
        ql=sqrt(1d0+(gammal+1d0)/2d0/gammal*(pstar/wl(5)-1d0))
        sl_temp=wl(2)-al*ql
    else
        !left rarefaction
        al=sqrt(gammal*wl(5)/wl(1))
        sl_temp=wl(2)-al
    end if
    if (pstar>wr(5)) then
        !right shock
        ar=sqrt(gammar*wr(5)/wr(1))
        !ar=sqrt(min(gammal,gammar)*wr(5)/wr(1))
        qr=sqrt(1d0+(gammar+1d0)/2d0/gammar*(pstar/wr(5)-1d0))
        sr_temp=wr(2)+ar*qr
    else
        !right rarefaction
        ar=sqrt(gammar*wr(5)/wr(1))
        sr_temp=wr(2)+ar
    end if
    sl=min(sl_temp,sr_temp)
    sr=max(sl_temp,sr_temp)
    sstar=(wr(5)-wl(5)+wl(1)*wl(2)*(sl-wl(2))-wr(1)*wr(2)*(sr-wr(2)))/(wl(1)*(sl-wl(2))-wr(1)*(sr-wr(2)))
end subroutine eos_hllc_analytic_slsrsstar

!subroutine eos_hllc_analytic_slsrsstar(wl,wr,sl,sr,sstar,egv,gamm,gammal,gammar,pstar)
!    !linearized estimate in "Restoration of the contact surface in the HLL-Riemann solver"
!    real(8), dimension(5) :: wl,wr
!    real(8) :: sl,sr,sl_temp,sr_temp,sstar,pstar,wsl,wsr,al,ar,ql,qr,gamm,gammal,gammar
!    real(8) :: sstar_real(2),p_jump,temp(2),egv(2),um,tml,tmr,aml,amr,rhol,rhor,rhot,at
!    rhot=half*(wl(1)+wr(1))
!    al=sqrt(gammal*wl(5)/wl(1))
!    ar=sqrt(gammar*wr(5)/wr(1))
!    at=half*(al+ar)
!    pstar=half*(wl(5)+wr(5))-half*(wr(2)-wl(2))*rhot*at
!    um=half*(wl(2)+wr(2))-half*(wr(5)-wl(5))/rhot/at
!    rhol=wl(1)+(wl(2)-um)*rhot/at
!    rhor=wr(1)+(um-wr(2))*rhot/at
!    print *,pstar,rhol,rhor
!    tml=solvetp(pstar,rhol)
!    tmr=solvetp(pstar,rhor)
!    print *,tml,tmr
!    aml=adiabatic_cs(rhol,tml)
!    amr=adiabatic_cs(rhor,tmr)
!    sl=min(wl(2)-al,um-aml)
!    sr=max(wr(2)+ar,um+amr)
!    !sstar=um
!    sstar=(wr(5)-wl(5)+wl(1)*wl(2)*(sl-wl(2))-wr(1)*wr(2)*(sr-wr(2)))/(wl(1)*(sl-wl(2))-wr(1)*(sr-wr(2)))
!end subroutine eos_hllc_analytic_slsrsstar

subroutine eos_hllc_analytic_ustar(w,egv,sstar,s,ustar)
    !given w, egv, sstar and sl or sr, calculate ustar
    real(8), dimension(5) :: w,ustar
    real(8) :: rho,t,sstar,h,egv,s,E,vx,vy,vz
    !egv is energy per unit volume (include radiation if present)
    !E is total energy
    h=w(1)*(s-w(2))/(s-sstar)
    ustar(1)=h
    ustar(2)=h*sstar
    ustar(3)=h*w(3)
    ustar(4)=h*w(4)
    rho=w(1)
    vx=w(2)
    vy=w(3)
    vz=w(4)
    E=egv+half*rho*(vx*vx+vy*vy+vz*vz)
    ustar(5)=h*(E/w(1)+(sstar-w(2))*(sstar+w(5)/w(1)/(s-w(2))))
end subroutine eos_hllc_analytic_ustar

subroutine eos_hllc_analytic_local_gamma(wl,wr,gamm,gammal,gammar,temp)
    !calculate gammas
    real(8), dimension(5) :: wl,wr
    real(8) :: gammal,gammar,p,rho,t
    real(8) :: gamm,pstar,e,rho_post,e_post,wsl,wsr,temp(2)
    integer :: i
    character(len=128) :: alert
    rho=wl(1)
    t=temp(1)
    call gammarhot(rho,t,gammal)
    rho=wr(1)
    t=temp(2)
    call gammarhot(rho,t,gammar)
    gamm=max(gammal,gammar)
end subroutine eos_hllc_analytic_local_gamma

subroutine eos_hllc_analytic_flux(wl,wr,temp,egv,flux,vlocalmax)
    !extended infomation:rho,vx,vy,vz,p,temp,egv
    !given left and right states (wl,wr), temp, and egv, calculate flux and vlocalmax
    real(8), dimension(5) :: wl,wr,flux,flux_temp,ur,ul,ustarl,ustarr
    real(8) :: sl,sr,sstar,vlocalmax,gamm,gammal,gammar,pstar,temp(2),egv(2)
    call eos_hllc_analytic_local_gamma(wl,wr,gamm,gammal,gammar,temp)
    call eos_hllc_analytic_slsrsstar(wl,wr,sl,sr,sstar,egv,gamm,gammal,gammar,pstar)
    vlocalmax=max(abs(sl),abs(sr))
    if (sl>0) then
        call eos_hllc_analytic_wtoxflux(wl,egv(1),flux)
    else if (sl<=0.and.sstar>0) then
        call eos_hllc_analytic_wtoxflux(wl,egv(1),flux_temp)
        call eos_hllc_analytic_wtou(wl,egv(1),ul)
        call eos_hllc_analytic_ustar(wl,egv(1),sstar,sl,ustarl)
        flux=flux_temp+sl*(ustarl-ul)
    else if (sstar<=0.and.sr>0) then
        call eos_hllc_analytic_wtoxflux(wr,egv(2),flux_temp)
        call eos_hllc_analytic_wtou(wr,egv(2),ur)
        call eos_hllc_analytic_ustar(wr,egv(2),sstar,sr,ustarr)
        flux=flux_temp+sr*(ustarr-ur)
    else
        call eos_hllc_analytic_wtoxflux(wr,egv(2),flux)
    end if
end subroutine eos_hllc_analytic_flux

subroutine eos_hllc_analytic_sweep(w,temp,egv,dudt,vglobalmax,dx,fluxarray,dir)
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
        call allocate_xsurface_data_block(vmaxarray)
        dudt=0d0
        do i=blk_xlb,blk_xub-1
            wl=w(i,1,1,:)
            wr=w(i+1,1,1,:)
            temp_local=(/temp(i,1,1),temp(i+1,1,1)/)
            egv_local=(/egv(i,1,1),egv(i+1,1,1)/)
            !print *,i,temp_local,wl(1),wl(5),wr(1),wr(5)
            call eos_hllc_analytic_flux(wl,wr,temp_local,egv_local,flux,vlocalmax)
            fluxarray(i,1,1,:)=flux
            vmaxarray(i,1,1)=vlocalmax
        end do
        vglobalmax=maxval(vmaxarray)
        do i=blk_xlb+1,blk_xub-1
            dudt(i,1,1,1:5)=(fluxarray(i-1,1,1,1:5)-fluxarray(i,1,1,1:5))/dx
        end do
        deallocate(vmaxarray)
    else if (nd==2) then
        !for 2d, the default is in x and y direction
        if (dir=='x') then
            call allocate_xsurface_data_block(vmaxarray)
            dudt=0d0
            do j=blk_ylb,blk_yub
                do i=blk_xlb,blk_xub-1
                    wl=w(i,j,1,:)
                    wr=w(i+1,j,1,:)
                    temp_local=(/temp(i,j,1),temp(i+1,j,1)/)
                    egv_local=(/egv(i,j,1),egv(i+1,j,1)/)
                    call eos_hllc_analytic_flux(wl,wr,temp_local,egv_local,flux,vlocalmax)
                    fluxarray(i,j,1,:)=flux
                    vmaxarray(i,j,1)=vlocalmax
                end do
            end do
            vglobalmax(1)=maxval(vmaxarray)
            do j=blk_ylb,blk_yub
                do i=blk_xlb+1,blk_xub-1
                    dudt(i,j,1,1:5)=(fluxarray(i-1,j,1,1:5)-fluxarray(i,j,1,1:5))/dx
                end do
            end do
            deallocate(vmaxarray)
        else if (dir=='y') then
            call allocate_ysurface_data(vmaxarray)
            dudt=0d0
            do j=blk_ylb,blk_yub-1
                do i=blk_xlb,blk_xub
                    wl=w(i,j,1,:)
                    wr=w(i,j+1,1,:)
                    temp_local=(/temp(i,j,1),temp(i,j+1,1)/)
                    egv_local=(/egv(i,j,1),egv(i,j+1,1)/)
                    call rotate_xyz_to_yzx(wl)
                    call rotate_xyz_to_yzx(wr)
                    call eos_hllc_analytic_flux(wl,wr,temp_local,egv_local,flux,vlocalmax)
                    call rotate_yzx_to_xyz(flux)
                    fluxarray(i,j,1,:)=flux
                    vmaxarray(i,j,1)=vlocalmax
                end do
            end do
            vglobalmax(2)=maxval(vmaxarray)
            do j=blk_ylb+1,blk_yub-1
                do i=blk_xlb,blk_xub
                    dudt(i,j,1,1:5)=(fluxarray(i,j-1,1,1:5)-fluxarray(i,j,1,1:5))/dx
                end do
            end do
            deallocate(vmaxarray)
        else
            alert='hydro_sweep no such direction'
            call abort_achilles(alert)
        end if
    else
    end if
end subroutine eos_hllc_analytic_sweep

subroutine hllc_muscl(blk)
    !use W_xl, W_xr, W_yl, W_yr to calculate fluxes
    type(blockdef), pointer :: blk
    real(8), dimension(5) :: wl,wr,flux
    real(8), dimension(:,:,:), allocatable :: vmaxarray
    real(8) :: vlocalmax,temp(2),egv(2),vblockmax(nd),dt(nd)
    integer :: i,j
    if (nd==1) then
        !for 1d the default is x direction
        blk%xflux=0d0
        call allocate_xsurface_data(vmaxarray)
        do i=0,blk_size_nx
            wl=blk%w_xr(i,1,1,1:5)
            wr=blk%w_xl(i+1,1,1,1:5)
            temp(1)=blk%temp_xr(i,1,1)
            temp(2)=blk%temp_xl(i+1,1,1)
            egv(1)=blk%egv_xr(i,1,1)
            egv(2)=blk%egv_xl(i+1,1,1)
            call eos_hllc_analytic_flux(wl,wr,temp,egv,flux,vlocalmax)
            blk%xflux(i,1,1,1:5)=flux
            vmaxarray(i,1,1)=vlocalmax
        end do
        if (time_sys%ntimestep==0) then
            vblockmax=maxval(vmaxarray)
            blk%dt_hydro_blk=blk%dxyz(1)*CFL/vblockmax(1)
        end if
        deallocate(vmaxarray)
    else if (nd==2) then
        !x direction first, then y direction
        dt=0d0
    end if
end subroutine hllc_muscl

subroutine eos_hllc_analytic_sweep_muscl(wl,wr,templ,tempr,egvl,egvr,vglobalmax,fluxarray,dir)
    !given the reconstructed wl, wr, templ, tempr, egvl, egvr, and dir, calculate fluxarray and vglobalmax
    real(8), dimension(:,:,:,:), allocatable :: wl,wr,fluxarray
    real(8), dimension(:,:,:), allocatable :: templ,tempr,egvl,egvr,vmaxarray
    real(8), dimension(5) :: w_left,w_right,flux
    real(8) :: vlocalmax,vglobalmax(nd),dx,temp_local(2),egv_local(2)
    character(len=1), optional :: dir
    character(len=128) :: alert
    integer :: i,j
    if (nd==1) then
        !for 1d, the default is x direction
        call allocate_xsurface_data_block(vmaxarray)
        do i=blk_xlb,blk_xub-1
            w_left=wr(i,1,1,1:5)
            w_right=wl(i+1,1,1,1:5)
            temp_local=(/tempr(i,1,1),templ(i+1,1,1)/)
            egv_local=(/egvr(i,1,1),egvl(i+1,1,1)/)
            call eos_hllc_analytic_flux(w_left,w_right,temp_local,egv_local,flux,vlocalmax)
            fluxarray(i,1,1,:)=flux
            vmaxarray(i,1,1)=vlocalmax
        end do
        vglobalmax=maxval(vmaxarray)
        deallocate(vmaxarray)
    else if (nd==2) then
        !for 2d, the default is in x and y direction
        if (dir=='x') then
            call allocate_xsurface_data_block(vmaxarray)
            do j=blk_ylb,blk_yub
                do i=blk_xlb,blk_xub-1
                    w_left=wr(i,j,1,1:5)
                    w_right=wl(i+1,j,1,1:5)
                    temp_local=(/tempr(i,j,1),templ(i+1,j,1)/)
                    egv_local=(/egvr(i,j,1),egvl(i+1,j,1)/)
                    call eos_hllc_analytic_flux(w_left,w_right,temp_local,egv_local,flux,vlocalmax)
                    fluxarray(i,j,1,1:5)=flux
                    vmaxarray(i,j,1)=vlocalmax
                end do
            end do
            vglobalmax(1)=maxval(vmaxarray)
            deallocate(vmaxarray)
        else if (dir=='y') then
            call allocate_ysurface_data_block(vmaxarray)
            do j=blk_ylb,blk_yub-1
                do i=blk_xlb,blk_xub
                    w_left=wr(i,j,1,1:5)
                    w_right=wl(i,j+1,1,1:5)
                    temp_local=(/tempr(i,j,1),templ(i,j+1,1)/)
                    egv_local=(/egvr(i,j,1),egvl(i,j+1,1)/)
                    call rotate_xyz_to_yzx(w_left)
                    call rotate_xyz_to_yzx(w_right)
                    call eos_hllc_analytic_flux(w_left,w_right,temp_local,egv_local,flux,vlocalmax)
                    call rotate_yzx_to_xyz(flux)
                    fluxarray(i,j,1,1:5)=flux
                    vmaxarray(i,j,1)=vlocalmax
                end do
            end do
            vglobalmax(2)=maxval(vmaxarray)
            deallocate(vmaxarray)
        else
            alert='hydro_sweep_muscl no such direction'
            call abort_achilles(alert)
        end if
    else
    end if
end subroutine eos_hllc_analytic_sweep_muscl

end module eos_hllc_analytic
