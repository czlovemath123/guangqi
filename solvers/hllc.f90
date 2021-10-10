module hllc
use phylib
use mathlib
use constant_gamma_eos
use datastructure
implicit none
contains

!function gleft(wl,p)
!    real(8), dimension(5) :: wl
!    real(8) :: p,al,bl,gleft
!    al=2/(gamma_gas+1)/wl(1)
!    bl=(gamma_gas-1)/(gamma_gas+1)*wl(5)
!    gleft=sqrt(al/(p+bl))
!end function gleft
!
!function gright(wr,p)
!    real(8), dimension(5) :: wr
!    real(8) :: p,ar,br,gright
!    ar=2/(gamma_gas+1)/wr(1)
!    br=(gamma_gas-1)/(gamma_gas+1)*wr(5)
!    gright=sqrt(ar/(p+br))
!end function gright
!
!subroutine hllc_pvrs(WL,WR,pstar)
!    real(8), dimension(5) :: WL,WR
!    real(8) :: ppvrs,rhoa,aa,al,ar,pstar
!    al=sqrt(gamma_gas*wl(5)/wl(1))
!    ar=sqrt(gamma_gas*wr(5)/wr(1))
!    aa=0.5*(al+ar)
!    rhoa=0.5*(WL(1)+WR(1))
!    ppvrs=0.5*(WL(5)+WR(5))-0.5*(WR(2)-WL(2))*rhoa*aa
!    pstar=max(0d0,ppvrs)
!end subroutine hllc_pvrs
!
!subroutine hllc_trrs(wl,wr,p)
!    real(8), dimension(5) :: wl,wr
!    real(8) :: p,z,al,ar
!    al=sqrt(gamma_gas*wl(5)/wl(1))
!    ar=sqrt(gamma_gas*wr(5)/wr(1))
!    z=(gamma_gas-1)/2/gamma_gas
!    p=pow((al+ar-gamma_gas*z*(wr(2)-wl(2)))/(al/pow(wl(5),z)+ar/pow(wr(5),z)),1/z)
!end subroutine hllc_trrs
!
!subroutine hllc_tsrs(wl,wr,p)
!    real(8), dimension(5) :: wl,wr
!    real(8) :: p,gl,gr,p0,pvrs
!    call hllc_pvrs(wl,wr,pvrs)
!    p0=max(real(0,kind=8),pvrs)
!    gl=gleft(wl,p0)
!    gr=gright(wr,p0)
!    p=(gl*wl(5)+gr*wr(5)-(wr(2)-wl(2)))/(gl+gr)
!end subroutine hllc_tsrs
!
!subroutine hllc_pstar(wl,wr,pstar)
!    real(8), dimension(5) :: wl,wr
!    real(8) :: pstar,q,qmax,pvrs,trrs,tsrs,pmax,pmin
!    call hllc_pvrs(wl,wr,pvrs)
!    pmax=max(wl(5),wr(5))
!    pmin=min(wl(5),wr(5))
!    q=pmax/pmin
!    qmax=2.0
!    if (q<=qmax.and.pvrs<=pmax.and.pvrs>=pmin) then
!        pstar=pvrs
!    else
!        if (pvrs<pmin) then
!            call hllc_trrs(wl,wr,pstar)
!        else
!            call hllc_tsrs(wl,wr,pstar)
!        end if
!    end if
!end subroutine hllc_pstar
!
!subroutine wavespeedestimate(WL,WR,pstar,Sstar,sl,sr)
!    real(8), dimension(5) :: WL,WR
!    real(8) :: sl_temp,sr_temp,sl,sr,ql,qr,pstar,Sstar,al,ar
!    al=sqrt(gamma_gas*WL(5)/WL(1))
!    ar=sqrt(gamma_gas*WR(5)/WR(1))
!    if (pstar<=WL(5)) then
!        ql=1
!    else
!        ql=sqrt(1+((1+gamma_gas)/(2*gamma_gas))*(pstar/WL(5)-1))
!    end if
!    if (pstar<=WR(5)) then
!        qr=1
!    else
!        qr=sqrt(1+((1+gamma_gas)/(2*gamma_gas))*(pstar/WR(5)-1))
!    end if
!    sl_temp=WL(2)-al*ql
!    sr_temp=WR(2)+ar*qr
!    sl=min(sl_temp,sr_temp)
!    sr=max(sl_temp,sr_temp)
!    Sstar=(WR(5)-WL(5)+WL(1)*WL(2)*(sl-WL(2))-WR(1)*WR(2)*(sr-WR(2)))/(WL(1)*(sl-WL(2))-WR(1)*(sr-WR(2)))
!    !print *,sl,sr,sstar,pstar
!end subroutine wavespeedestimate
!
!subroutine calUstar(W,Ustar,Sstar,s)
!    real(8), dimension(5) :: W,Ustar
!    real(8) :: Sstar,s,Energy,h,vx,vy,vz,v
!    h=W(1)*(s-W(2))/(s-Sstar)
!    ustar(1)=h
!    ustar(2)=h*Sstar
!    ustar(3)=h*w(3)
!    ustar(4)=h*w(4)
!    vx=w(2)
!    vy=w(3)
!    vz=w(4)
!    v=sqrt(vx*vx+vy*vy+vz*vz)
!    Energy=W(1)*(1/(gamma_gas-1)*W(5)/W(1)+0.5*(vx*vx+vy*vy+vz*vz))
!    Ustar(5)=h*(Energy/W(1)+(Sstar-W(2))*(Sstar+W(5)/(W(1)*(s-W(2)))))
!end subroutine calUstar
!
!subroutine hllc_flux(wl,wr,flux,vlocalmax)
!    !given left and right states (wl,wr), calculate flux and vlocalmax
!    real(8), dimension(5) :: WL,WR,flux,Ustar,UL,UR,fl,fr
!    real(8) :: Sstar,sl,sr,pstar,vlocalmax
!    fl=0;fr=0
!    call hllc_pstar(WL,WR,pstar)
!    call wavespeedestimate(WL,WR,pstar,Sstar,sl,sr)
!    vlocalmax=max(abs(sl),abs(sr))
!    if (sl>=0) then
!        call wtoflux(WL,flux)
!    else if (sl<0 .and. Sstar>=0) then
!        call WtoU(WL,UL)
!        call wtoflux(WL,fl)
!        call calUstar(WL,Ustar,Sstar,sl)
!        flux=sl*(Ustar-UL)+fl
!    else if (Sstar<0 .and. sr>=0) then
!        call wtou(WR,UR)
!        call wtoflux(WR,fr)
!        call calUstar(WR,Ustar,Sstar,sr)
!        flux=sr*(Ustar-UR)+fr
!    else
!        call wtoflux(WR,flux)
!    end if
!    !print *,flux
!end subroutine hllc_flux


!P. Batten, N. Clarke, C. Lambert, and D. M. Causon "On the choince of wavespeeds for the HLLC Riemann solver" 1997

subroutine hllc_roe_averaged(wl,wr,q_roe,c_roe)
    !use Roe average Riemann solver to find q_roe, c_roe
    real(8), dimension(5) :: wl,wr,w_roe
    real(8) :: egv(2),hl,hr,rhot,ut,vt,wt,ht,h,tt,q_roe,c_roe
    real(8) :: sqrtrhol,sqrtrhor,sqrtrho,ep
    egv=(/wl(5)/(gamma_gas-1d0),wr(5)/(gamma_gas-1d0)/)
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
    ep=(ht-half*(ut**2+vt**2+wt**2))
    c_roe=sqrt((gamma_gas-1d0)*ep)
    q_roe=ut
end subroutine hllc_roe_averaged

function hllc_p(wl,wr)
    real(8) :: wl(5),wr(5),hllc_p
    real(8) :: sl,sr,sm,vlocalmax,egv(2)
    real(8) :: ql,qr,q_roe,cl,cr,c_roe,dsl,dsr,dsm
    real(8) :: rhom,pm,momum,momvm,momwm,el,er,em
    call hllc_roe_averaged(wl,wr,q_roe,c_roe)
    ql=wl(2)
    qr=wr(2)
    cl=sqrt(gamma_gas*wl(5)/wl(1))
    cr=sqrt(gamma_gas*wr(5)/wr(1))
    egv(1)=wl(5)/(gamma_gas-1d0)
    egv(2)=wr(5)/(gamma_gas-1d0)
    sl=min(ql-cl,q_roe-c_roe)
    sr=max(qr+cr,q_roe+c_roe)
    dsl=sl-ql
    dsr=sr-qr
    sm=(wr(1)*qr*dsr-wl(1)*ql*dsl+wl(5)-wr(5))/(wr(1)*dsr-wl(1)*dsl)
    if (sl>zero) then
        hllc_p=wl(5)
    else if (sl<=zero.and.sm>zero) then
        dsm=sl-sm
        rhom=wl(1)*dsl/dsm
        pm=wl(5)-wl(1)*dsl*(ql-sm)
        momum=(dsl*wl(1)*wl(2)+(pm-wl(5)))/dsm
        momvm=dsl*wl(1)*wl(3)/dsm
        momwm=dsl*wl(1)*wl(4)/dsm
        el=egv(1)+half*wl(1)*(wl(2)**2+wl(3)**2+wl(4)**2)
        em=(dsl*el-ql*wl(5)+pm*sm)/dsm
        hllc_p=pm
    else if (sm<=zero.and.sr>zero) then
        dsm=sr-sm
        rhom=wr(1)*dsr/dsm
        pm=wr(5)-wr(1)*dsr*(qr-sm)
        momum=(dsr*wr(1)*wr(2)+(pm-wr(5)))/dsm
        momvm=dsr*wr(1)*wr(3)/dsm
        momwm=dsr*wr(1)*wr(4)/dsm
        er=egv(2)+half*wr(1)*(wr(2)**2+wr(3)**2+wr(4)**2)
        em=(dsr*er-qr*wr(5)+pm*sm)/dsm
        hllc_p=pm
    else
        hllc_p=wr(5)
    end if
end function hllc_p

subroutine hllc_flux(wl,wr,flux,vlocalmax,debug)
    !given left and right states (wl,wr), calculate flux and vlocalmax
    real(8), dimension(5) :: wl,wr,flux,ul,ur,fl,fr
    real(8) :: sl,sr,sm,vlocalmax,egv(2)
    real(8) :: ql,qr,q_roe,cl,cr,c_roe,dsl,dsr,dsm
    real(8) :: rhom,pm,momum,momvm,momwm,el,er,em
    logical, optional :: debug
    call hllc_roe_averaged(wl,wr,q_roe,c_roe)
    ql=wl(2)
    qr=wr(2)
    cl=sqrt(gamma_gas*wl(5)/wl(1))
    cr=sqrt(gamma_gas*wr(5)/wr(1))
    egv(1)=wl(5)/(gamma_gas-1d0)
    egv(2)=wr(5)/(gamma_gas-1d0)
    sl=min(ql-cl,q_roe-c_roe)
    sr=max(qr+cr,q_roe+c_roe)
    dsl=sl-ql
    dsr=sr-qr
    sm=(wr(1)*qr*dsr-wl(1)*ql*dsl+wl(5)-wr(5))/(wr(1)*dsr-wl(1)*dsl)
    if (present(debug)) then
        print *,'haha',sm,wr(1)*qr*dsr,wl(1)*ql*dsl,wl(5),wr(5),(wr(1)*qr*dsr-wl(1)*ql*dsl+wl(5)-wr(5))
    end if
    if (sl>zero) then
        call wtoflux(wl,flux)
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
    else if (sm==0) then
        call wtou(wl,ul)
        call wtou(wr,ur)
        call wtoflux(wl,fl)
        call wtoflux(wr,fr)
        flux=(sr*fl-sl*fr+sl*sr*(ur-ul))/(sr-sl)
    else if (sm<zero.and.sr>=zero) then
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
        call wtoflux(wr,flux)
    end if
    vlocalmax=max(abs(sl),abs(sr))
end subroutine hllc_flux

subroutine hllc_sweep(w,vblockmax,dx,fluxarray,dir,reverse)
    !given w and dir, calculate dudt and vblockmax
    !u is the conserved quantities in dudt
    real(8), dimension(:,:,:,:), allocatable :: w
    real(8), allocatable :: vmaxarray(:,:,:),fluxarray(:,:,:,:)
    real(8), dimension(5) :: wl,wr,flux
    real(8) :: vlocalmax,vblockmax(nd),dx
    character(len=1), optional :: dir
    logical, optional :: reverse
    character(len=128) :: alert
    integer :: i,j
    !dx=dy=dz is the computational size of the control volume
    if (nd==1) then
        !for 1d, the default is x direction
        call allocate_xsurface_data_block(vmaxarray)
        do i=blk_xlb,blk_xub-1
            wl=w(i,1,1,:)
            wr=w(i+1,1,1,:)
            call hllc_flux(wl,wr,flux,vlocalmax)
            fluxarray(i,1,1,:)=flux
            vmaxarray(i,1,1)=vlocalmax
        end do
        vblockmax=maxval(vmaxarray)
        deallocate(vmaxarray)
    else if (nd==2) then
        !for 2d, the default is x and y direction
        if (dir=='x') then
            call allocate_xsurface_data_block(vmaxarray)
            do j=blk_ylb,blk_yub
                do i=blk_xlb,blk_xub-1
                    wl=w(i,j,1,:)
                    wr=w(i+1,j,1,:)
                    call hllc_flux(wl,wr,flux,vlocalmax)
                    fluxarray(i,j,1,:)=flux
                    vmaxarray(i,j,1)=vlocalmax
                end do
            end do
            vblockmax(1)=maxval(vmaxarray)
            deallocate(vmaxarray)
        else if (dir=='y') then
            call allocate_ysurface_data_block(vmaxarray)
            do j=blk_ylb,blk_yub-1
                do i=blk_xlb,blk_xub
                    wl=w(i,j,1,:)
                    wr=w(i,j+1,1,:)
                    if (igeometry==2) then
                        if (reverse) then
                            wl(3)=-wl(3)
                            wr(3)=-wr(3)
                            call rotate_xyz_to_yzx(wl)
                            call rotate_xyz_to_yzx(wr)
                            call hllc_flux(wr,wl,flux,vlocalmax)
                            call rotate_yzx_to_xyz(flux)
                            flux=-flux
                            flux(3)=-flux(3)
                            fluxarray(i,j,1,1:5)=flux
                        else
                            call rotate_xyz_to_yzx(wl)
                            call rotate_xyz_to_yzx(wr)
                            call hllc_flux(wl,wr,flux,vlocalmax)
                            call rotate_yzx_to_xyz(flux)
                            fluxarray(i,j,1,1:5)=flux
                        end if
                    else
                        call rotate_xyz_to_yzx(wl)
                        call rotate_xyz_to_yzx(wr)
                        call hllc_flux(wl,wr,flux,vlocalmax)
                        call rotate_yzx_to_xyz(flux)
                        fluxarray(i,j,1,1:5)=flux
                    end if
                    vmaxarray(i,j,1)=vlocalmax
                end do
            end do
            vblockmax(2)=maxval(vmaxarray)
            deallocate(vmaxarray)
        end if
    else
    end if
end subroutine hllc_sweep

subroutine hllc_muscl(blk)
    !use W_xl, W_xr, W_yl, W_yr to calculate fluxes
    type(blockdef), pointer :: blk
    real(8), dimension(5) :: wl,wr,flux
    real(8), dimension(:,:,:), allocatable :: vmaxarray
    real(8) :: vlocalmax,vblockmax,dt(nd)
    integer :: i,j
    if (nd==1) then
        !for 1d, the default is x direction
        blk%xflux=0d0
        call allocate_xsurface_data(vmaxarray)
        do i=0,blk_size_nx
            wl=blk%w_xr(i,1,1,1:5)
            wr=blk%w_xl(i+1,1,1,1:5)
            call hllc_flux(wl,wr,flux,vlocalmax)
            blk%xflux(i,1,1,1:5)=flux
            vmaxarray(i,1,1)=vlocalmax
        end do
        if (time_sys%ntimestep==0) then
            vblockmax=maxval(vmaxarray)
            blk%dt_hydro_blk=blk%dxyz(1)*CFL/vblockmax
        end if
        deallocate(vmaxarray)
    else if (nd==2) then
        !x direction first, then y direction
        dt=0d0;blk%xflux=0d0;blk%yflux=0d0
        call allocate_xsurface_data_block(vmaxarray)
        do j=1,blk_size_ny
            do i=0,blk_size_nx
                wl=blk%w_xr(i,j,1,1:5)
                wr=blk%w_xl(i+1,j,1,1:5)
                call hllc_flux(wl,wr,flux,vlocalmax)
                blk%xflux(i,j,1,1:5)=flux
                vmaxarray(i,j,1)=vlocalmax
            end do
        end do
        if (time_sys%ntimestep==0) then
            vblockmax=maxval(vmaxarray)
            dt(1)=blk%dxyz(1)*CFL/vblockmax
        end if
        deallocate(vmaxarray)
        call allocate_ysurface_data_block(vmaxarray)
        do j=0,blk_size_ny
            do i=1,blk_size_nx
                wl=blk%w_yr(i,j,1,1:5)
                wr=blk%w_yl(i,j+1,1,1:5)
                call rotate_xyz_to_yzx(wl)
                call rotate_xyz_to_yzx(wr)
                call hllc_flux(wl,wr,flux,vlocalmax)
                call rotate_yzx_to_xyz(flux)
                blk%yflux(i,j,1,1:5)=flux
                vmaxarray(i,j,1)=vlocalmax
            end do
        end do
        if (time_sys%ntimestep==0) then
            vblockmax=maxval(vmaxarray)
            dt(2)=blk%dxyz(2)*CFL/vblockmax
            blk%dt_hydro_blk=minval(dt)
        end if
        deallocate(vmaxarray)
    end if
end subroutine hllc_muscl

end module hllc
