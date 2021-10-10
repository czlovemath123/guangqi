module exactriemann
use mathlib
use phylib
implicit none
contains
!below: calculate pstar
function extrmfl(x)
    real(8), dimension(:), allocatable, intent(in) :: x
    real(8) :: AL,BL,cs,extrmfl
    !x(1) will be pstar, x(2:5) is WL
    AL=2/((x(5)+1)*x(2))
    BL=(x(5)-1)/(x(5)+1)*x(4)
    cs=sqrt(x(5)*x(4)/x(2))
    if (x(1)>x(4)) then
        extrmfl=(x(1)-x(4))*sqrt(AL/(x(1)+BL))
    else
        extrmfl=2*cs/(x(5)-1)*((x(1)/x(4))**((x(5)-1)/(2*x(5)))-1)
    end if
end function extrmfl

function extrmfr(x)
    real(8), dimension(:), allocatable, intent(in) :: x
    real(8) :: AR,BR,cs,extrmfr
    !x(1) will be pstar, x(2:5) is WR
    AR=2/((x(5)+1)*x(2))
    BR=(x(5)-1)/(x(5)+1)*x(4)
    cs=sqrt(x(5)*x(4)/x(2))
    if (x(1)>x(4)) then
        extrmfr=(x(1)-x(4))*sqrt(AR/(x(1)+BR))
    else
        extrmfr=2*cs/(x(5)-1)*((x(1)/x(4))**((x(5)-1)/(2*x(5)))-1)
    end if
end function extrmfr

function extrmflderi(x)
    real(8), dimension(:), allocatable, intent(in) :: x
    real(8) :: AL,BL,cs,extrmflderi
    !x(1) will be pstar, x(2:5) is WL
    AL=2/((x(5)+1)*x(2))
    BL=(x(5)-1)/(x(5)+1)*x(4)
    cs=sqrt(x(5)*x(4)/x(2))
    if (x(1)>x(4)) then
        extrmflderi=sqrt(AL/(BL+x(1)))*(1-(x(1)-x(4))/(2*(BL+x(1))))
    else
        extrmflderi=1/(x(2)*cs)*((x(1)/x(4))**(-(x(5)+1)/(2*x(5))))
    end if
end function extrmflderi

function extrmfrderi(x)
    real(8), dimension(:), allocatable, intent(in) :: x
    real(8) :: AR,BR,cs,extrmfrderi
    !x(1) will be pstar, x(2:5) is WR
    AR=2/((x(5)+1)*x(2))
    BR=(x(5)-1)/(x(5)+1)*x(4)
    cs=sqrt(x(5)*x(4)/x(2))
    if (x(1)>x(4)) then
        extrmfrderi=sqrt(AR/(BR+x(1)))*(1-(x(1)-x(4))/(2*(BR+x(1))))
    else
        extrmfrderi=1/(x(2)*cs)*(x(1)/x(4))**(-(x(5)+1)/(2*x(5)))
    end if
end function extrmfrderi

function extrmf(x)
    real(8),  dimension(:), allocatable :: x,xl,xr
    !x(1) will be pstar, x(2:5) is WL, x(6:9) is WR
    real(8) :: extrmf
    allocate(xl(5),xr(5))
    xl(1)=x(1)
    xl(2:5)=x(2:5)
    xr(1)=x(1)
    xr(2:5)=x(6:9)
    extrmf=extrmfl(xl)+extrmfr(xr)+(xr(3)-xl(3))
    deallocate(xl,xr)
end function extrmf

function extrmfderi(x)
    real(8),  dimension(:), allocatable :: x,xl,xr
    !x(1) will be pstar, x(2:5) is WL, x(6:9) is WR
    real(8) :: extrmfderi
    allocate(xl(5),xr(5))
    xl(1)=x(1)
    xl(2:5)=x(2:5)
    xr(1)=x(1)
    xr(2:5)=x(6:9)
    extrmfderi=extrmflderi(xl)+extrmfrderi(xr)
    deallocate(xl,xr)
end function extrmfderi

!given pstar, calculate ustarl and ustarr, to get the p-u phase curve
function exactrm_ustarl(pstar,wl)
    real(8) :: wl(4),pstar,exactrm_ustarl
    real(8), allocatable :: x(:)
    allocate(x(5))
    x(1)=pstar
    x(2:5)=wl
    exactrm_ustarl=wl(2)-extrmfl(x)
    deallocate(x)
end function exactrm_ustarl

function exactrm_ustarr(pstar,wr)
    real(8) :: wr(4),pstar,exactrm_ustarr
    real(8), allocatable :: x(:)
    allocate(x(5))
    x(1)=pstar
    x(2:5)=wr
    exactrm_ustarr=wr(2)+extrmfr(x)
    deallocate(x)
end function exactrm_ustarr

!below: calculate 5 nonlinear wave velocities

function extrmustar(wl,wr,pstar)
    real(8) :: wl(4),wr(4),pstar,extrmustar
    real(8), dimension(:), allocatable :: xl,xr
    allocate(xl(7),xr(7))
    xl(1)=pstar
    xl(2:5)=wl
    xr(1)=pstar
    xr(2:5)=wr
    extrmustar=0.5*(wl(2)+wr(2))+0.5*(extrmfr(xr)-extrmfl(xl))
    deallocate(xl,xr)
end function extrmustar

function extrmvleftshock(wl,pstar)
    real(8) :: wl(4),pstar,extrmvleftshock,csl
    csl=sqrt(wl(4)*wl(3)/wl(1))
    extrmvleftshock=wl(2)-csl*sqrt((wl(4)+1)/(2*wl(4))*pstar/wl(3)+(wl(4)-1)/(2*wl(4)))
end function extrmvleftshock

subroutine extrmvleftrare(wl,pstar,ustar,shl,stl)
    real(8) :: wl(4),pstar,ustar,shl,stl,csl,csstar
    csl=sqrt(wl(4)*wl(3)/wl(1))
    csstar=csl*(pstar/wl(3))**((wl(4)-1)/(2*wl(4)))
    shl=wl(2)-csl
    stl=ustar-csstar
end subroutine extrmvleftrare

function extrmvrightshock(wr,pstar)
    real(8) :: wr(4),pstar,extrmvrightshock,csr
    csr=sqrt(wr(4)*wr(3)/wr(1))
    extrmvrightshock=wr(2)+csr*sqrt((wr(4)+1)/(2*wr(4))*pstar/wr(3)+(wr(4)-1)/(2*wr(4)))
end function extrmvrightshock

subroutine extrmvrightrare(wr,pstar,ustar,shl,stl)
    real(8) :: wr(4),pstar,ustar,shl,stl,csr,csstar
    csr=sqrt(wr(4)*wr(3)/wr(1))
    csstar=csr*(pstar/wr(3))**((wr(4)-1)/(2*wr(4)))
    shl=wr(2)+csr
    stl=ustar+csstar
end subroutine extrmvrightrare

!****************************************************************
!below: calculate 4 nonlinear states

subroutine extrmleftshock(wl,w,pstar,ustar)
    real(8), dimension(4) :: wl,w
    real(8) :: pstar,rhostar,ustar
    rhostar=wl(1)*(pstar/wl(3)+(wl(4)-1)/(wl(4)+1))/((wl(4)-1)/(wl(4)+1)*pstar/wl(3)+1)
    w(1)=rhostar
    w(2)=ustar
    w(3)=pstar
    w(4)=wl(4)
end subroutine extrmleftshock

subroutine extrmleftrarestar(wl,w,pstar,ustar)
    real(8), dimension(4) :: wl,w
    real(8) :: pstar,rhostar,ustar
    rhostar=wl(1)*(pstar/wl(3))**(1/wl(4))
    w(1)=rhostar
    w(2)=ustar
    w(3)=pstar
    w(4)=wl(4)
end subroutine extrmleftrarestar

subroutine extrmleftrarefan(wl,w,pstar,x,t)
    real(8), dimension(4) :: wl,w
    real(8) :: pstar,x,rho,u,p,t,csl
    !x is the coordinate, x/t is the value we need. t is in info
    csl=sqrt(wl(4)*wl(3)/wl(1))
    rho=wl(1)*(2/(wl(4)+1)+(wl(4)-1)/((wl(4)+1)*csl)*(wl(2)-x/t))**(2/(wl(4)-1))
    u=2/(wl(4)+1)*(csl+(wl(4)-1)/2*wl(2)+x/t)
    p=wl(3)*(2/(wl(4)+1)+(wl(4)-1)/((wl(4)+1)*csl)*(wl(2)-x/t))**(2*wl(4)/(wl(4)-1))
    w(1)=rho
    w(2)=u
    w(3)=p
    w(4)=wl(4)
end subroutine extrmleftrarefan

subroutine extrmrightshock(wr,w,pstar,ustar)
    real(8), dimension(4) :: wr,w
    real(8) :: pstar,rhostar,ustar
    rhostar=wr(1)*(pstar/wr(3)+(wr(4)-1)/(wr(4)+1))/((wr(4)-1)/(wr(4)+1)*pstar/wr(3)+1)
    w(1)=rhostar
    w(2)=ustar
    w(3)=pstar
    w(4)=wr(4)
end subroutine extrmrightshock

subroutine extrmrightrarestar(wr,w,pstar,ustar)
    real(8), dimension(4) :: wr,w
    real(8) :: pstar,rhostar,ustar
    rhostar=wr(1)*(pstar/wr(3))**(1/wr(4))
    w(1)=rhostar
    w(2)=ustar
    w(3)=pstar
    w(4)=wr(4)
end subroutine extrmrightrarestar

subroutine extrmrightrarefan(wr,w,pstar,x,t)
    real(8), dimension(4) :: wr,w
    real(8) :: pstar,x,rho,u,p,t,csr
    !x is the coordinate, x/t is the value we need. t is in info
    csr=sqrt(wr(4)*wr(3)/wr(1))
    rho=wr(1)*(2/(wr(4)+1)-(wr(4)-1)/((wr(4)+1)*csr)*(wr(2)-x/t))**(2/(wr(4)-1))
    u=2/(wr(4)+1)*(-csr+(wr(4)-1)/2*wr(2)+x/t)
    p=wr(3)*(2/(wr(4)+1)-(wr(4)-1)/((wr(4)+1)*csr)*(wr(2)-x/t))**(2*wr(4)/(wr(4)-1))
    w(1)=rho
    w(2)=u
    w(3)=p
    w(4)=wr(4)
end subroutine extrmrightrarefan

function extrmg(AL,BL,p)
    real(8) :: AL,BL,p,extrmg
    extrmg=sqrt(AL/(p+BL))
end function extrmg


!******************************************************************
!below: sample the solutions

subroutine extrmsample(wl,wr,x,w,tolerance,t)
    real(8), dimension(4) :: wl,wr,w
    real(8) :: pstar,ustar,x,pguess,tolerance,csl,csr,p0,AL,BL,AR,BR,t
    real(8) :: vleftshock,vrightshock,vleftshl,vleftstl,vrightshl,vrightstl
    real(8), dimension(:), allocatable :: estimate
    procedure(fun), pointer :: ptr1,ptr2
    allocate(estimate(13))
    csl=sqrt(wl(4)*wl(3)/wl(1))
    csr=sqrt(wr(4)*wr(3)/wr(1))
    if ((wr(2)-wl(2)).ge.(csl+csr)) then
        pguess=((csl+csr-0.5*(wl(4)-1)*(wr(2)-wl(2)))/(csl/(wl(3)**((wl(4)-1)/(2*wl(4)))) &
            +csr/(wr(3)**((wr(4)-1)/(2*wr(4))))))**(2*wl(4)/(wl(4)-1))
    else
        AL=2/(wl(4)+1)/wl(1)
        BL=(wl(4)-1)/(wl(4)+1)*wl(3)
        AR=2/(wr(4)+1)/wr(1)
        BR=(wr(4)-1)/(wr(4)+1)*wr(3)
        p0=0.5*(wl(3)+wr(3))-0.125*(wr(2)-wl(2))*(wr(1)+wl(1))*(csr+csl)
        p0=max(tolerance,p0)
        pguess=(extrmg(AL,BL,p0)*wl(3)+extrmg(AR,BR,p0)*wr(3)-(wr(2)-wl(2)))/(extrmg(AL,BL,p0)+extrmg(AR,BR,p0))
        pguess=max(tolerance,pguess)
    end if
    estimate(1)=pguess
    estimate(2:5)=wl
    estimate(6:9)=wr
    ptr1=>extrmf
    ptr2=>extrmfderi
    call newtonraphson(ptr1,ptr2,pstar,estimate,tolerance)
    w(1)=pstar
    !print *,'****************This is exact Riemann solver******************************'
    !print *,'wl= ',wl
    !print *,'wr= ',wr
    !print *,'pstar= ',pstar
    !print *,'****************This is exact Riemann solver******************************'
    deallocate(estimate)
    ustar=extrmustar(wl,wr,pstar)
    if (x/t<ustar) then              !on the left side
        if (pstar>wl(3)) then        !left shock
            vleftshock=extrmvleftshock(wl,pstar)
            if (x/t<vleftshock) then
                w=wl
            else
                call extrmleftshock(wl,w,pstar,ustar)
            end if
        else                         !left rarefaction
            call extrmvleftrare(wl,pstar,ustar,vleftshl,vleftstl)
            if (x/t<vleftshl) then
                w=wl
            else if (x/t>vleftstl) then
                call extrmleftrarestar(wl,w,pstar,ustar)
            else
                call extrmleftrarefan(wl,w,pstar,x,t)
            end if
        end if
    else                             !on the right side
        if (pstar>wr(3)) then        !right shock
            vrightshock=extrmvrightshock(wr,pstar)
            if (x/t>vrightshock) then
                w=wr
            else
                call extrmrightshock(wr,w,pstar,ustar)
            end if
        else                         !left rarefaction
            call extrmvrightrare(wr,pstar,ustar,vrightshl,vrightstl)
            if (x/t>vrightshl) then
                w=wr
            else if (x/t<vrightstl) then
                call extrmrightrarestar(wr,w,pstar,ustar)
            else
                call extrmrightrarefan(wr,w,pstar,x,t)
            end if
        end if
    end if
end subroutine extrmsample

end module exactriemann
