module eos_hllc_tabulated
use datastructure
use eos_tabulated
use mathlib
use phylib
implicit none
contains
subroutine eos_hllc_pvrs(wl,wr,pvrs,gammal,gammar)
    real(8), dimension(5) :: wl,wr
    real(8) :: pvrs,al,ar,gammal,gammar
    al=sqrt(gammal*wl(5)/wl(1))  !eos_tabulated_sound_speed(wl(1),wl(5))
    ar=sqrt(gammar*wr(5)/wr(1))  !eos_tabulated_sound_speed(wr(1),wr(5))
    pvrs=half*(wl(5)+wr(5))-half*(wr(2)-wl(2))*half*(wl(1)+wr(1))*half*(al+ar)
end subroutine eos_hllc_pvrs

subroutine eos_hllc_trrs(wl,wr,p,gammal,gammar,gamm)
    real(8), dimension(5) :: wl,wr
    real(8) :: p,z,al,ar,gammal,gammar,gamm
    z=(gamm-1)/two/gamm
    al=sqrt(gammal*wl(5)/wl(1))  !eos_tabulated_sound_speed(wl(1),wl(5))
    ar=sqrt(gammar*wr(5)/wr(1))  !eos_tabulated_sound_speed(wr(1),wr(5))
    p=pow((al+ar-gamm*z*(wr(2)-wl(2)))/(al/pow(wl(5),z)+ar/pow(wr(5),z)),1/z)
end subroutine eos_hllc_trrs

function eos_gleft(wl,p,gamm)
    real(8), dimension(5) :: wl
    real(8) :: p,al,bl,eos_gleft,gamm
    al=2/(gamm+1)/wl(1)
    bl=(gamm-1)/(gamm+1)*wl(5)
    eos_gleft=sqrt(al/(p+bl))
end function eos_gleft

function eos_gright(wr,p,gamm)
    real(8), dimension(5) :: wr
    real(8) :: p,ar,br,eos_gright,gamm
    ar=2/(gamm+1)/wr(1)
    br=(gamm-1)/(gamm+1)*wr(5)
    eos_gright=sqrt(ar/(p+br))
end function eos_gright

subroutine eos_hllc_tsrs(wl,wr,p,gammal,gammar)
    real(8), dimension(5) :: wl,wr
    real(8) :: p,gl,gr,p0,pvrs,gammal,gammar
    call eos_hllc_pvrs(wl,wr,pvrs,gammal,gammar)
    p0=max(zero,pvrs)
    gl=eos_gleft(wl,p0,gammal)
    gr=eos_gright(wr,p0,gammar)
    p=(gl*wl(5)+gr*wr(5)-(wr(2)-wl(2)))/(gl+gr)
end subroutine eos_hllc_tsrs

subroutine eos_hllc_pstar(wl,wr,pstar,gammal,gammar,gamm)
    real(8), dimension(5) :: wl,wr
    real(8) :: pstar,q,qmax,pvrs,trrs,tsrs,pmax,pmin,gammal,gammar,gamm
    call eos_hllc_pvrs(wl,wr,pvrs,gammal,gammar)
    pmax=max(wl(5),wr(5))
    pmin=min(wl(5),wr(5))
    q=pmax/pmin
    qmax=2.0
    if (q<=qmax.and.pvrs<=pmax.and.pvrs>=pmin) then
        pstar=pvrs
    else
        if (pvrs<pmin) then
            call eos_hllc_trrs(wl,wr,pstar,gammal,gammar,gamm)
        else
            call eos_hllc_tsrs(wl,wr,pstar,gammal,gammar)
        end if
    end if
end subroutine eos_hllc_pstar

subroutine eos_hllc_slsrsstar(wl,wr,sl,sr,sstar,gammal,gammar,gamm,pstar)
    real(8), dimension(5) :: wl,wr,w
    real(8) :: sl,sr,sl_temp,sr_temp,sstar,pstar,pstar2,wsl,wsr,al,ar,ql,qr,gamm
    real(8) :: sstar_real(2),gammal,gammar
    procedure(fun), pointer :: ptr
    integer :: i
    call eos_hllc_pstar(wl,wr,pstar,gammal,gammar,gamm)
    !calculate sl,sr and sstar
    if (pstar>wl(5)) then
        !left shock
        al=sqrt(gammal*wl(5)/wl(1))
        ql=sqrt(1+(gammal+1)/2/gammal*(pstar/wl(5)-1))
        sl_temp=wl(2)-al*ql
    else
        !left rarefaction
        al=sqrt(gammal*wl(5)/wl(1))
        sl_temp=wl(2)-al
    end if
    if (pstar>wr(5)) then
        !right shock
        ar=sqrt(gammar*wr(5)/wr(1))
        qr=sqrt(1+(gammar+1)/2/gammar*(pstar/wr(5)-1))
        sr_temp=wr(2)+ar*qr
    else
        !right rarefaction
        ar=sqrt(gammar*wr(5)/wr(1))
        sr_temp=wr(2)+ar
    end if
    sl=min(sl_temp,sr_temp)
    sr=max(sl_temp,sr_temp)
    sstar=(wr(5)-wl(5)+wl(1)*wl(2)*(sl-wl(2))-wr(1)*wr(2)*(sr-wr(2)))/(wl(1)*(sl-wl(2))-wr(1)*(sr-wr(2)))
end subroutine eos_hllc_slsrsstar

subroutine eos_hllc_ustar(w,ustar,sstar,s,pstar)
    real(8), dimension(5) :: w,ustar,w_dxdt
    real(8) :: rho,sstar,pstar,h,eg,s,E,vx,vy,vz,v,dxdt
    !eg is energy per unit mass
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
    v=sqrt(vx*vx+vy*vy+vz*vz)
    eg=eos_tabulated_internal_e_per_mass(w(1),w(5))
    E=rho*eg+half*rho*v*v
    ustar(5)=h*(E/w(1)+(sstar-w(2))*(sstar+w(5)/w(1)/(s-w(2))))
end subroutine eos_hllc_ustar

subroutine local_gamma(wl,wr,gammal,gammar,gamm)
    real(8), dimension(5) :: wl,wr
    real(8) :: gammal,gammar,p,rho,temp,log10_pressure,log10_rho,log10_temp,ktoe
    real(8) :: gamm,gamm_new,pstar,du,el,er,elmax,ermax
    real(8) :: e,rho_post,e_post,wsl,wsr
    procedure(fun), pointer :: ptr
    integer :: i
    !estimate the pstar
    rho=wl(1)
    p=wl(5)
    gammal=eos_tabulated_gamma(rho,p)
    rho=wr(1)
    p=wr(5)
    gammar=eos_tabulated_gamma(rho,p)
    gamm=max(gammal,gammar)
end subroutine local_gamma

subroutine eos_hllc_flux(wl,wr,flux,vlocalmax)
    real(8), dimension(5) :: wl,wr,w
    real(8), dimension(5) :: flux,flux_temp,ur,ul,ustarl,ustarr
    real(8) :: sl,sr,sstar,pstar,vlocalmax,gamm,gammal,gammar
    call local_gamma(wl,wr,gammal,gammar,gamm)
    call eos_hllc_slsrsstar(wl,wr,sl,sr,sstar,gammal,gammar,gamm,pstar)
    vlocalmax=max(abs(sl),abs(sr))
    if (sl>0) then
        call eos_wtoxflux(wl,flux)
    else if (sl<=0.and.sstar>0) then
        call eos_wtoxflux(wl,flux_temp)
        call eos_wtou(wl,ul)
        call eos_hllc_ustar(wl,ustarl,sstar,sl,pstar)
        flux=flux_temp+sl*(ustarl-ul)
    else if (sstar<=0.and.sr>0) then
        call eos_wtoxflux(wr,flux_temp)
        call eos_wtou(wr,ur)
        call eos_hllc_ustar(wr,ustarr,sstar,sr,pstar)
        flux=flux_temp+sr*(ustarr-ur)
    else
        call eos_wtoxflux(wr,flux)
    end if
end subroutine eos_hllc_flux

subroutine eos_hllc_sweep(w,dudt,vglobalmax,dx)
    real(8), dimension(:,:,:,:), allocatable, intent(in) :: w
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: dudt
    real(8), allocatable :: vmaxarray(:,:,:),fluxarray(:,:,:,:)
    real(8), dimension(5) :: wl,wr,flux
    real(8) :: vlocalmax,vglobalmax(nd),dx
    character(len=128) :: alert
    integer :: i
    allocate(vmaxarray(nx+1,1,1),fluxarray(nx+1,1,1,5))
    vmaxarray=0d0
    fluxarray=0d0
    dudt=0d0
    do i=0,nx
        wl=w(i,1,1,:)
        wr=w(i+1,1,1,:)
        call eos_hllc_flux(wl,wr,flux,vlocalmax)
        fluxarray(i+1,1,1,:)=flux
        vmaxarray(i+1,1,1)=vlocalmax
    end do
    vglobalmax=maxval(vmaxarray)
    do i=1,nx
        dudt(i,1,1,1:5)=(fluxarray(i,1,1,1:5)-fluxarray(i+1,1,1,1:5))/dx
    end do
    deallocate(vmaxarray,fluxarray)
end subroutine eos_hllc_sweep

end module eos_hllc_tabulated
