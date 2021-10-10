module eosriemannsolver
use mathlib
use phylib
use eos_tabulated
use exactriemann
use datastructure
implicit none
contains

function ws2_root_fun(x)
    real(8), dimension(:), allocatable :: x
    !x(1) ws2, x(2) rho, x(3) pressure, x(4) pstar, x(5) gamm
    real(8) :: ws2_root_fun,rhostar,gamm,pstar,estar,e,ratio
    pstar=x(4)
    gamm=x(5)
    rhostar=1/(1/x(2)-(pstar-x(3))/x(1))
    ratio=x(4)/x(3)
    !print *,'rhostar=',rhostar
    if (rhostar.lt.rho_min_bound) rhostar=min(rho_max_bound,ratio*x(2))
    if (rhostar.gt.rho_max_bound) rhostar=rho_max_bound
    !print *,'x=',x,rhostar
    estar=eos_tabulated_internal_e_per_mass(rhostar,pstar)
    e=eos_tabulated_internal_e_per_mass(x(2),x(3))
    ws2_root_fun=x(1)*(estar-e)-0.5*(pstar**2-x(3)**2)
    !print *,ws2_root_fun,x(1),rhostar
end function ws2_root_fun

subroutine eos_exact_wsl_ext(wl,wr,ws2,pstar,gamm)
    real(8), dimension(4) :: wl,wr
    real(8) :: pstar,ws2,rho,p,rho_post,gamm
    rho=wl(1)
    p=wl(3)
    !call eos_riemann_pstar_ext(wl,wr,pstar,gamm)
    rho_post=rho*(pstar/p+(gamm-1)/(gamm+1))/((gamm-1)/(gamm+1)*pstar/p+1)
    ws2=-(pstar-p)/(1/rho_post-1/rho)
end subroutine eos_exact_wsl_ext

subroutine eos_exact_wsr_ext(wl,wr,ws2,pstar,gamm)
    real(8), dimension(4) :: wl,wr
    real(8) :: pstar,ws2,rho,p,rho_post,gamm
    rho=wr(1)
    p=wr(3)
    !call eos_riemann_pstar_ext(wl,wr,pstar,gamm)
    rho_post=rho*(pstar/p+(gamm-1)/(gamm+1))/((gamm-1)/(gamm+1)*pstar/p+1)
    ws2=-(pstar-p)/(1/rho_post-1/rho)
end subroutine eos_exact_wsr_ext

subroutine find_xleft_ws_speed(ptr,wl,wr,pstar,ws)
    real(8), dimension(4) :: wl,wr
    real(8), dimension(:), allocatable :: x,root
    real(8) :: pstar,ws,ws2min,ws2max,ws2,csl,csr,converge
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    !ws is the mass flux across shock
    !ptr points to defined ws2_root_fun function
    conv_mode=2
    converge=2d0*epsilon(wl(1))
    csr=eos_tabulated_sound_speed(wr(1),wr(3))
    csl=eos_tabulated_sound_speed(wl(1),wl(3))
    !call eos_exact_wsl_ext(wl,wr,ws2max,1.667d0)
    !call eos_exact_wsl_ext(wl,wr,ws2min,1.01d0)
    !ws2max=ws2max*10
    !ws2min=ws2min/10
    ws2min=0d0
    ws2max=(wr(1)*(abs(wr(2))+csr)*100)**2
    allocate(x(5),root(5))
    x(1)=ws2max/2d0
    x(2)=wl(1)
    x(3)=wl(3)
    x(4)=pstar
    x(5)=wl(4)
    if (pstar>wl(3)*1.0001) then
        call bisectionroot(ptr,ws2min,ws2max,x,converge,root,conv_mode)
    else
        call eos_exact_wsl_ext(wl,wr,root(1),pstar,1.4d0)
    end if
    ws=sqrt(root(1))
    deallocate(x,root)
end subroutine find_xleft_ws_speed

subroutine find_xright_ws_speed(ptr,wl,wr,pstar,ws)
    real(8), dimension(4) :: wl,wr
    real(8), dimension(:), allocatable :: x,root
    real(8) :: pstar,ws,ws2min,ws2max,ws2,csl,csr,converge
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    !ws is the mass flux across shock
    !ptr points to defined find_ws_root function
    conv_mode=2
    converge=2d0*epsilon(wl(1))
    csl=eos_tabulated_sound_speed(wl(1),wl(3))
    csr=eos_tabulated_sound_speed(wr(1),wr(3))
    !call eos_exact_wsr_ext(wl,wr,ws2max,1.667d0)
    !call eos_exact_wsr_ext(wl,wr,ws2min,1.01d0)
    !ws2max=ws2max*10
    !ws2min=ws2min/10
    ws2min=0d0
    ws2max=(wl(1)*(abs(wl(2))+csl)*100)**2
    allocate(x(5),root(5))
    x(1)=ws2max/2d0
    x(2)=wr(1)
    x(3)=wr(3)
    x(4)=pstar
    x(5)=wr(4)
    if (pstar>wr(3)*1.0001) then
        call bisectionroot(ptr,ws2min,ws2max,x,converge,root,conv_mode)
    else
        call eos_exact_wsr_ext(wl,wr,root(1),pstar,1.4d0)
    end if
    ws=sqrt(root(1))
    deallocate(x,root)
end subroutine find_xright_ws_speed

subroutine find_xleft_rare_star(pstar,ustar,rhostar,wl)
    real(8) :: pstar,ustar,pressure,rhonew,rhoold,c,dp,unew,uold,rhostar
    real(8), dimension(4) :: wl
    integer :: i,n
    real(8) :: log10_pressure,log10_temp,log10_rho,gamm,dp_finer,p_step,fp
    !c is Lagrangian wave speed
    pressure=wl(3)
    rhoold=wl(1)
    uold=wl(2)
    n=max(ceiling((log10(pressure)-log10(pstar))*20000),1000)
    !print *,n,pressure/pstar
    !integration point
    dp=(pstar-pressure)/n
    !dp_finer=dp/20
    c=rhoold*eos_tabulated_sound_speed(rhoold,pressure)
    rhonew=1d0/(1d0/rhoold-1d0/c/c*dp)
    unew=uold-1d0/c*dp
    pressure=pressure+dp
    do i=1,n
        c=rhoold*eos_tabulated_sound_speed(rhoold,pressure)
        rhonew=1d0/(1d0/rhoold-1d0/c/c*dp)
        unew=uold-1d0/c*dp
        rhoold=rhonew
        uold=unew
        pressure=pressure+dp
    end do
    ustar=unew
    rhostar=rhonew
end subroutine find_xleft_rare_star

subroutine find_xright_rare_star(pstar,ustar,rhostar,wr)
    real(8) :: pstar,ustar,pressure,rhonew,rhoold,c,dp,unew,uold,rhostar
    real(8), dimension(4) :: wr
    integer :: i,n
    real(8) :: log10_pressure,log10_temp,log10_rho,gamm,dp_finer,p_step,fp
    !c is Lagrangian wave speed
    pressure=wr(3)
    rhoold=wr(1)
    uold=wr(2)
    n=max(ceiling((log10(pressure)-log10(pstar))*20000),1000)
    !integration point
    dp=(pstar-pressure)/n
    !dp_finer=dp/20
    c=rhoold*eos_tabulated_sound_speed(rhoold,pressure)
    rhonew=1d0/(1d0/rhoold-1d0/c/c*dp)
    unew=uold+1d0/c*dp
    pressure=pressure+dp
    do i=1,n
        c=rhoold*eos_tabulated_sound_speed(rhoold,pressure)
        rhonew=1d0/(1d0/rhoold-1d0/c/c*dp)
        unew=uold+1d0/c*dp
        rhoold=rhonew
        uold=unew
        pressure=pressure+dp
    end do
    ustar=unew
    rhostar=rhonew
end subroutine find_xright_rare_star

subroutine find_xleft_rare_fan(pstar,wl,dxdt,w)
    real(8) :: pstar,pressure,rhonew,rhoold,c,csold,csnew,dp,unew,uold,rhostar,dxdt
    real(8) :: du,fp,rhofan,ufan,pfan,log10_pressure,log10_temp,log10_rho,gamm,dp_finer,p_step
    real(8), dimension(4) :: wl,w
    integer :: i,n
    !c is Lagrangian wave speed
    n=5e5   !max(ceiling((log10(pressure)-log10(pstar))*20000),1000)
    rhoold=wl(1)
    uold=wl(2)
    pressure=wl(3)
    csold=eos_tabulated_sound_speed(rhoold,pressure)
    c=rhoold*csold
    dp=(pstar-pressure)/n
    !dp_finer=dp/20
    rhonew=1/(1/rhoold-1/c/c*dp)
    unew=uold-1/c*dp
    pressure=pressure+dp
    csnew=eos_tabulated_sound_speed(rhonew,pressure)
    do while (unew-csnew.lt.dxdt)
        rhoold=rhonew
        uold=unew
        csold=csnew
        c=rhoold*csold   !eos_tabulated_sound_speed(rhoold,pressure)
        p_step=dp
        rhonew=1d0/(1d0/rhoold-1d0/c/c*p_step)
        unew=uold-1d0/c*p_step
        pressure=pressure+p_step
        csnew=eos_tabulated_sound_speed(rhonew,pressure)
    end do
    du=dxdt-uold+csold
    fp=du/(-1d0/c)
    rhofan=1d0/(1d0/rhoold-1d0/c/c*fp)
    ufan=uold-1d0/c*fp
    pfan=pressure-p_step+fp
    w=(/rhofan,ufan,pfan,wl(4)/)
end subroutine find_xleft_rare_fan

subroutine find_xright_rare_fan(pstar,wr,dxdt,w)
    real(8) :: pstar,pressure,rhonew,rhoold,c,csold,csnew,dp,unew,uold,rhostar,dxdt
    real(8) :: du,fp,rhofan,ufan,pfan,p_step
    real(8), dimension(4) :: wr,w
    integer :: n
    !c is the Lagrangian wave speed
    n=5e5   !max(ceiling((log10(pressure)-log10(pstar))*20000),1000)
    rhoold=wr(1)
    uold=wr(2)
    pressure=wr(3)
    csold=eos_tabulated_sound_speed(rhoold,pressure)
    c=rhoold*csold
    dp=(pstar-pressure)/n
    !dp_finer=dp/20
    rhonew=1/(1/rhoold-1/c/c*dp)
    unew=uold+1/c*dp
    pressure=pressure+dp
    csnew=eos_tabulated_sound_speed(rhonew,pressure)
    do while (unew+csnew.gt.dxdt)
        rhoold=rhonew
        uold=unew
        csold=csnew
        c=rhoold*csold   !eos_tabulated_sound_speed(rhoold,pressure)
        p_step=dp
        rhonew=1/(1/rhoold-1/c/c*p_step)
        unew=uold+1/c*p_step
        pressure=pressure+p_step
        csnew=eos_tabulated_sound_speed(rhonew,pressure)
    end do
    du=dxdt-uold-csold
    fp=du/(1/c)
    rhofan=1/(1/rhoold-1/c/c*fp)
    ufan=uold+1/c*fp
    pfan=pressure-p_step+fp
    w=(/rhofan,ufan,pfan,wr(4)/)
end subroutine find_xright_rare_fan

function eos_exactrm_ustarl(pstar,wl,wr)
    !given a pstar and wl, (wr is only for informational purpose), calculate ustarl with bisection method
    real(8) :: pstar,eos_exactrm_ustarl,wl(4),wr(4)
    real(8), dimension(:), allocatable :: x,root
    real(8) :: rhostar,ws,ws2min,ws2max,ws2,csl,csr,converge,ws2max1,ws2max2
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    if (pstar>wl(3)) then
        !ws is the mass flux across shock
        !ptr points to defined ws2_root_fun function
        ptr=>ws2_root_fun
        conv_mode=2
        converge=2d0*epsilon(wl(1))
        csl=eos_tabulated_sound_speed(wl(1),wl(3))
        csr=eos_tabulated_sound_speed(wr(1),wr(3))
        ws2min=0d0
        ws2max1=(wl(1)*(abs(wl(2))+csl)*100)**2
        ws2max2=(wr(1)*(abs(wr(2))+csr)*100)**2
        ws2max=max(ws2max1,ws2max2)
        allocate(x(5),root(5))
        x(1)=ws2max/2d0
        x(2)=wl(1)
        x(3)=wl(3)
        x(4)=pstar
        x(5)=wl(4)
        if (pstar>wl(3)*1.0001) then
            call bisectionroot(ptr,ws2min,ws2max,x,converge,root,conv_mode)
        else
            call eos_exact_wsl_ext(wl,wr,root(1),pstar,1.4d0)
        end if
        ws=sqrt(root(1))
        deallocate(x,root)
        eos_exactrm_ustarl=wl(2)-(pstar-wl(3))/ws
    else
        call find_xleft_rare_star(pstar,eos_exactrm_ustarl,rhostar,wl)
    end if
end function eos_exactrm_ustarl

function eos_exactrm_ustarr(pstar,wl,wr)
    !given a pstar and wr, (wl is only for informational purpose), calculate ustarr with bisection method
    real(8) :: pstar,eos_exactrm_ustarr,wl(4),wr(4)
    real(8), dimension(:), allocatable :: x,root
    real(8) :: rhostar,ws,ws2min,ws2max,ws2,csl,csr,converge,ws2max1,ws2max2
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    if (pstar>wr(3)) then
        !ws is the mass flux across shock
        !ptr points to defined find_ws_root function
        ptr=>ws2_root_fun
        conv_mode=2
        converge=2d0*epsilon(wl(1))
        csl=eos_tabulated_sound_speed(wl(1),wl(3))
        csr=eos_tabulated_sound_speed(wr(1),wr(3))
        ws2min=0d0
        ws2max1=(wl(1)*(abs(wl(2))+csl)*100)**2
        ws2max2=(wr(1)*(abs(wr(2))+csr)*100)**2
        ws2max=max(ws2max1,ws2max2)
        allocate(x(5),root(5))
        x(1)=ws2max/2d0
        x(2)=wr(1)
        x(3)=wr(3)
        x(4)=pstar
        x(5)=wr(4)
        if (pstar>wr(3)*1.0001) then
            call bisectionroot(ptr,ws2min,ws2max,x,converge,root,conv_mode)
        else
            call eos_exact_wsr_ext(wl,wr,root(1),pstar,1.4d0)
        end if
        ws=sqrt(root(1))
        deallocate(x,root)
        eos_exactrm_ustarr=wr(2)+(pstar-wr(3))/ws
    else
        call find_xright_rare_star(pstar,eos_exactrm_ustarr,rhostar,wr)
    end if
end function eos_exactrm_ustarr

subroutine eos_riemann_pstar_ext(wl,wr,pstar,gamm)
    real(8), dimension(4) :: wl,wr
    real(8) :: pstar,gamm,csl,csr,pguess,AL,BL,AR,BR,p0,tolerance
    real(8), dimension(:), allocatable :: estimate
    procedure(fun), pointer :: ptr1,ptr2
    allocate(estimate(9))
    tolerance=1e-12
    csl=sqrt(gamm*wl(3)/wl(1))
    csr=sqrt(gamm*wr(3)/wr(1))
    if ((wr(2)-wl(2)).ge.(csl+csr)) then
        pguess=((csl+csr-0.5*(gamm-1)*(wr(2)-wl(2)))/(csl/(wl(3)**((gamm-1)/(2*gamm))) &
            +csr/(wr(3)**((gamm-1)/(2*gamm)))))**(2*gamm/(gamm-1))
    else
        AL=2/(gamm+1)/wl(1)
        BL=(gamm-1)/(gamm+1)*wl(3)
        AR=2/(gamm+1)/wr(1)
        BR=(gamm-1)/(gamm+1)*wr(3)
        p0=0.5*(wl(3)+wr(3))-0.125*(wr(2)-wl(2))*(wr(1)+wl(1))*(csr+csl)
        p0=max(tolerance,p0)
        pguess=(extrmg(AL,BL,p0)*wl(3)+extrmg(AR,BR,p0)*wr(3)-(wr(2)-wl(2)))/(extrmg(AL,BL,p0)+extrmg(AR,BR,p0))
        pguess=max(tolerance,pguess)
    end if
    estimate(1)=pguess
    estimate(2:4)=wl(1:3)
    estimate(5)=gamm
    estimate(6:8)=wr(1:3)
    estimate(9)=gamm
    ptr1=>extrmf
    ptr2=>extrmfderi
    call newtonraphson(ptr1,ptr2,pstar,estimate,tolerance)
    if (isnan(pstar)) pstar=10**p_min_bound
    deallocate(estimate)
end subroutine eos_riemann_pstar_ext

subroutine secant_wsl(p,ul,wl,wr,wsl)
    real(8) :: p,ul,wl(4),wr(4),wsl
    real(8), dimension(:), allocatable :: x,root
    real(8) :: rhostar,ws,ws2min,ws2max,ws2,csl,csr,converge,ws2max1,ws2max2
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    if (p>wl(3)) then
        ptr=>ws2_root_fun
        conv_mode=2
        converge=2d0*epsilon(wl(1))
        csl=eos_tabulated_sound_speed(wl(1),wl(3))
        csr=eos_tabulated_sound_speed(wr(1),wr(3))
        ws2min=0d0
        ws2max1=(wl(1)*(abs(wl(2))+csl)*100)**2
        ws2max2=(wr(1)*(abs(wr(2))+csr)*100)**2
        ws2max=max(ws2max1,ws2max2)
        allocate(x(5),root(5))
        x(1)=ws2max/2d0
        x(2)=wl(1)
        x(3)=wl(3)
        x(4)=p
        x(5)=wl(4)
        if (p>wl(3)*1.0001) then
            call bisectionroot(ptr,ws2min,ws2max,x,converge,root,conv_mode)
        else
            call eos_exact_wsl_ext(wl,wr,root(1),p,1.4d0)
        end if
        wsl=sqrt(root(1))
        deallocate(x,root)
    else
        if (ul==wl(2)) then
            wsl=eos_tabulated_sound_speed(wl(1),wl(3))
        else
            wsl=abs(p-wl(3))/abs(ul-wl(2))
        end if
    end if
end subroutine secant_wsl

subroutine secant_wsr(p,ur,wl,wr,wsr)
    real(8) :: p,ur,wl(4),wr(4),wsr
    real(8), dimension(:), allocatable :: x,root
    real(8) :: rhostar,ws,ws2min,ws2max,ws2,csl,csr,converge,ws2max1,ws2max2
    procedure(fun), pointer :: ptr
    integer :: conv_mode
    if (p>wr(3)) then
        ptr=>ws2_root_fun
        conv_mode=2
        converge=2d0*epsilon(wl(1))
        csl=eos_tabulated_sound_speed(wl(1),wl(3))
        csr=eos_tabulated_sound_speed(wr(1),wr(3))
        ws2min=0d0
        ws2max1=(wl(1)*(abs(wl(2))+csl)*100)**2
        ws2max2=(wr(1)*(abs(wr(2))+csr)*100)**2
        ws2max=max(ws2max1,ws2max2)
        allocate(x(5),root(5))
        x(1)=ws2max/2d0
        x(2)=wr(1)
        x(3)=wr(3)
        x(4)=p
        x(5)=wr(4)
        if (p>wr(3)*1.0001) then
            call bisectionroot(ptr,ws2min,ws2max,x,converge,root,conv_mode)
        else
            call eos_exact_wsr_ext(wl,wr,root(1),p,1.4d0)
        end if
        wsr=sqrt(root(1))
        deallocate(x,root)
    else
        if (ur==wr(2)) then
            wsr=eos_tabulated_sound_speed(wr(1),wr(3))
        else
            wsr=abs(p-wr(3))/abs(ur-wr(2))
        end if
    end if
end subroutine secant_wsr

function ml(wl,pstar)
    real(8) :: wl(4),ml,pstar,w
    w=pstar/wl(3)
    if (w>=1) then
        ml=sqrt(wl(1)*wl(3))*sqrt((wl(4)+1)/2*w+(wl(4)-1)/2)
    else
        ml=sqrt(wl(1)*wl(3))*((wl(4)-1)/2/sqrt(wl(4))*(1-w)/(1-w**((wl(4)-1)/2/wl(4))))
    end if
end function ml

function mr(wr,pstar)
    real(8) :: wr(4),mr,pstar,w
    w=pstar/wr(3)
    if (w>=1) then
        mr=sqrt(wr(1)*wr(3))*sqrt((wr(4)+1)/2*w+(wr(4)-1)/2)
    else
        mr=sqrt(wr(1)*wr(3))*((wr(4)-1)/2/sqrt(wr(4))*(1-w)/(1-w**((wr(4)-1)/2/wr(4))))
    end if
end function mr

subroutine eos_riemann_secant(wl,wr,pstar,ustar)
    real(8), dimension(4) :: wl,wr
    real(8) :: pstar,ustar,p0,p1,p2,ul0,ul1,ur0,ur1,wsl0,wsl1,wsr0,wsr1
    real(8), dimension(:), allocatable :: x,root,y
    real(8) :: rhostar,ws,ws2min,ws2max,ws2,csl,csr,converge,ws2max1,ws2max2
    real(8) :: tolerance,p_temp,AL,BL,AR,BR
    procedure(fun), pointer :: ptr
    integer :: conv_mode,i
    i=0

    !allocate(y(13))
    !tolerance=1e-10
    !csl=sqrt(wl(4)*wl(3)/wl(1))
    !csr=sqrt(wr(4)*wr(3)/wr(1))
    !if ((wr(2)-wl(2))>=(csl+csr)) then
    !    p0=((csl+csr-0.5*(wl(4)-1)*(wr(2)-wl(2)))/(csl/(wl(3)**((wl(4)-1)/(2*wl(4)))) &
    !        +csr/(wr(3)**((wr(4)-1)/(2*wr(4))))))**(2*wl(4)/(wl(4)-1))
    !else
    !    AL=2/(wl(4)+1)/wl(1)
    !    BL=(wl(4)-1)/(wl(4)+1)*wl(3)
    !    AR=2/(wr(4)+1)/wr(1)
    !    BR=(wr(4)-1)/(wr(4)+1)*wr(3)
    !    p_temp=0.5*(wl(3)+wr(3))-0.125*(wr(2)-wl(2))*(wr(1)+wl(1))*(csr+csl)
    !    p_temp=max(tolerance,p_temp)
    !    p0=(extrmg(AL,BL,p0)*wl(3)+extrmg(AR,BR,p0)*wr(3)-(wr(2)-wl(2)))/(extrmg(AL,BL,p0)+extrmg(AR,BR,p0))
    !    p0=max(tolerance,p0)
    !end if
    !y(1)=p0
    !y(2:7)=wl
    !y(8:13)=wr
    !p1=p0-extrmf(y)/extrmfderi(y)
    !deallocate(y)

    p0=(wl(3)+wr(3))/2
    p1=(wl(2)-wr(2)+wl(3)/ml(wl,p0)+wr(3)/mr(wr,p0))/(1d0/ml(wl,p0)+1d0/mr(wr,p0))
    p1=max(1d-10,p1)

    !call eos_riemann_pstar_ext(wl,wr,p0,1.667d0)
    !call eos_riemann_pstar_ext(wl,wr,p1,1.05d0)
    !p0=p0/2
    !p1=p1*2

    ul0=exactrm_ustarl(p0,wl)
    ul1=exactrm_ustarl(p1,wl)
    ur0=exactrm_ustarr(p0,wr)
    ur1=exactrm_ustarr(p1,wr)
    do while (abs(p1-p0)>1e-9)
        call secant_wsl(p0,ul0,wl,wr,wsl0)
        call secant_wsl(p1,ul1,wl,wr,wsl1)
        call secant_wsr(p0,ur0,wl,wr,wsr0)
        call secant_wsr(p1,ur1,wl,wr,wsr1)
        ul0=wl(2)-(p0-wl(3))/wsl0
        ur0=wr(2)+(p0-wr(3))/wsr0
        ul1=wl(2)-(p1-wl(3))/wsl1
        ur1=wr(2)+(p1-wr(3))/wsr1
        p2=p1-(ur1-ul1)*abs(p1-p0)/(abs(ul1-ul0)+abs(ur1-ur0))
        p0=p1
        p1=p2
        i=i+1
        print *,i,p0,ul1,ur1
    end do
    pstar=p1
    ustar=(ul1+ur1)/2
end subroutine eos_riemann_secant

subroutine eos_riemann(wl,wr,pstar,ustar,rhostar,lambda)
    real(8), dimension(4) :: wl,wr
    real(8) :: pstar,ustar
    real(8) :: ustarl(2),ustarr(2),wsl,wsr,pressure(2),csl,csr, &
    gamm,pressure_new,ustarl_new,ustarr_new,mach,pguess1,pguess2
    procedure(fun), pointer :: ptr,ptr1,ptr2
    real(8), dimension(:), allocatable :: estimate
    integer :: i,slope
    real(8), dimension(2) :: lambda,rhostar
    logical :: mark
    !lambda(1) is the fastest left wave speed
    !lambda(2) is the fastest right wave speed
    !rho(1) is the left density in star region
    !rho(2) is the right density in star region
    !print *,'********************This is EOS common pstar exact solver*************************'
    gamm=(wl(1)*wl(4)+wr(1)*wr(4))/(wl(1)+wr(1))
    allocate(estimate(13))
    csl=eos_tabulated_sound_speed(wl(1),wl(3))
    csr=eos_tabulated_sound_speed(wr(1),wr(3))
    if ((wr(2)-wl(2)).gt.0) then
        call eos_riemann_pstar_ext(wl,wr,pguess1,1.667d0)
        call eos_riemann_pstar_ext(wl,wr,pguess2,1.05d0)
        !pressure=(/pguess1,pguess2/)
        pressure(1)=max(10**p_min_bound,pguess1/2)
        pressure(2)=min(max(wl(3),wr(3)),pguess2*2)
        mark=.true.
        !the upper bound of pstar is no more than the current max pressure
    else
        call eos_riemann_pstar_ext(wl,wr,pguess1,1.667d0)
        call eos_riemann_pstar_ext(wl,wr,pguess2,1.05d0)
        !pressure=(/pguess1,pguess2/)
        pressure(1)=max(min(wl(3),wr(3)),pguess1/2)
        if (pguess2.gt.10**p_max_bound) then
            print *,'intermediate state pressure may exceed the max'
        end if
        pressure(2)=min(10**p_max_bound,pguess2*2)
        !print *,pressure
        mark=.false.
        !the lower bound of pstar is no less than the current min pressure
    end if
    do i=1,2
        if (pressure(i)>wl(3)) then
            ptr=>ws2_root_fun
            call find_xleft_ws_speed(ptr,wl,wr,pressure(i),wsl)
            ustarl(i)=wl(2)-(pressure(i)-wl(3))/wsl
            lambda(1)=wl(2)-wsl/wl(1)
            rhostar(1)=1d0/(1d0/wl(1)-(pressure(2)-wl(3))/wsl**2)
        else
            call find_xleft_rare_star(pressure(i),ustarl(i),rhostar(1),wl)
            lambda(1)=wl(2)-csl
        end if
        if (pressure(i)>wr(3)) then
            ptr=>ws2_root_fun
            call find_xright_ws_speed(ptr,wl,wr,pressure(i),wsr)
            ustarr(i)=wr(2)+(pressure(i)-wr(3))/wsr
            lambda(2)=wr(2)+wsr/wr(1)
            rhostar(2)=1d0/(1d0/wr(1)-(pressure(2)-wr(3))/wsr**2)
            !print *,wsr
        else
            call find_xright_rare_star(pressure(i),ustarr(i),rhostar(2),wr)
            lambda(2)=wr(2)+csr
        end if
    end do
    !print *,'estimate',ustarl
    !print *,'estimate',ustarr
    if (mark) then
        if ((ustarr(2)-ustarl(2))*(ustarr(1)-ustarl(1)).gt.0) then
            print *, 'bisection wrong, solution not bracketted, two rare'
            stop
        else if ((ustarr(2)-ustarl(2)).gt.(ustarr(1)-ustarl(1))) then
            slope=1
        else
            slope=-1
        end if
        do while((abs(pressure(2)-pressure(1)).gt.2d0*epsilon(abs(pressure(1))+abs(pressure(2)))))
            pressure_new=(pressure(1)+pressure(2))/2d0
            if (pressure_new>wl(3)) then
                ptr=>ws2_root_fun
                call find_xleft_ws_speed(ptr,wl,wr,pressure_new,wsl)
                ustarl_new=wl(2)-(pressure_new-wl(3))/wsl
                lambda(1)=wl(2)-wsl/wl(1)
                rhostar(1)=1d0/(1d0/wl(1)-(pressure_new-wl(3))/wsl**2)
                !print *,wsl,pressure_new,rhostar(1)
            else
                call find_xleft_rare_star(pressure_new,ustarl_new,rhostar(1),wl)
                lambda(1)=wl(2)-csl
            end if
            if (pressure_new>wr(3)) then
                ptr=>ws2_root_fun
                call find_xright_ws_speed(ptr,wl,wr,pressure_new,wsr)
                ustarr_new=wr(2)+(pressure_new-wr(3))/wsr
                lambda(2)=wr(2)+wsr/wr(1)
                rhostar(2)=1d0/(1d0/wr(1)-(pressure_new-wr(3))/wsr**2)
                !print *,wsr,pressure_new,rhostar(2)
            else
                call find_xright_rare_star(pressure_new,ustarr_new,rhostar(2),wr)
                lambda(2)=wr(2)+csr
            end if
            if ((ustarr_new-ustarl_new)*slope.gt.0) then
                pressure(2)=pressure_new
            else
                pressure(1)=pressure_new
            end if
            i=i+1
            print *,i,pressure_new,ustarl_new,ustarr_new
        !    print *,'**************************'
        !    print *,ustarl_new,ustarr_new,pressure_new/wl(3)
        end do
        pstar=(pressure(1)+pressure(2))/2
        ustar=(ustarl_new+ustarr_new)/2
        !print *,'ustar',ustarl_new,ustarr_new
        !print *,'pressure diff',abs(pressure(2)-pressure(1))
        !print *,'pstar',pstar
    else
        !print *,ustarr(2)-ustarl(2),ustarr(1)-ustarl(1)
        if ((ustarr(2)-ustarl(2))*(ustarr(1)-ustarl(1)).gt.0) then
            print *, 'bisection wrong, solution not bracketted, two shock'
            stop
        else if ((ustarr(2)-ustarl(2)).gt.(ustarr(1)-ustarl(1))) then
            slope=1
        else
            slope=-1
        end if
        do while(abs(pressure(2)-pressure(1)).gt.1e-9.or.(abs(ustarr_new-ustarl_new).gt.1e-3))
            pressure_new=(pressure(1)+pressure(2))/2d0
            if (pressure_new>wl(3)) then
                ptr=>ws2_root_fun
                call find_xleft_ws_speed(ptr,wl,wr,pressure_new,wsl)
                ustarl_new=wl(2)-(pressure_new-wl(3))/wsl
                lambda(1)=wl(2)-wsl/wl(1)
                rhostar(1)=1/(1/wl(1)-(pressure_new-wl(3))/wsl**2)
                !print *,wsl,pressure_new,rhostar(1)
            else
                call find_xleft_rare_star(pressure_new,ustarl_new,rhostar(1),wl)
                lambda(1)=wl(2)-csl
            end if
            if (pressure_new>wr(3)) then
                ptr=>ws2_root_fun
                call find_xright_ws_speed(ptr,wl,wr,pressure_new,wsr)
                ustarr_new=wr(2)+(pressure_new-wr(3))/wsr
                lambda(2)=wr(2)+wsr/wr(1)
                rhostar(2)=1/(1/wr(1)-(pressure_new-wr(3))/wsr**2)
                !print *,wsr,pressure_new,rhostar(2)
            else
                call find_xright_rare_star(pressure_new,ustarr_new,rhostar(2),wr)
                lambda(2)=wr(2)+csr
            end if
            if ((ustarr_new-ustarl_new)*slope.gt.0) then
                pressure(2)=pressure_new
            else
                pressure(1)=pressure_new
            end if
            i=i+1
            print *,i,pressure_new,ustarl_new,ustarr_new
            !print *,'**************************'
            !print *,ustarl_new,ustarr_new,pressure_new/wl(3)
        end do
        pstar=(pressure(1)+pressure(2))/2
        ustar=(ustarl_new+ustarr_new)/2
        !print *,'ustar',ustarl_new,ustarr_new
        !print *,'pressure diff',abs(pressure(2)-pressure(1))
        !print *,'pstar',pstar
    end if
    !print *,'lambda=',lambda
    !print *,'result',ustarl_new,ustarr_new
    !print *,'********************This is EOS common pstar exact solver*************************'
end subroutine eos_riemann

subroutine eos_xstate(wl,wr,w,dxdt,pstar,ustar,rhostar,lambda)
    real(8), dimension(4) :: wl,wr,w
    real(8) :: dxdt
    real(8) :: csstar
    real(8), optional :: pstar,ustar,rhostar(2),lambda(2)
    if (present(pstar)) then
        !print *,'Riemann problem data fed'
    else
        call eos_riemann(wl,wr,pstar,ustar,rhostar,lambda)
    end if
    !print *,pstar,ustar,lambda
    if (dxdt.lt.lambda(1)) then
        w=wl
    else if (dxdt.lt.ustar .and. dxdt.ge.lambda(1)) then
        if (pstar.gt.wl(3)) then    !left shock
            w(1)=rhostar(1)
            w(2)=ustar
            w(3)=pstar
            w(4)=wl(4)
        else                        !left rarefaction
            csstar=eos_tabulated_sound_speed(rhostar(1),pstar)
            if (dxdt.lt.ustar-csstar) then
                call find_xleft_rare_fan(pstar,wl,dxdt,w)
            else
                w(1)=rhostar(1)
                w(2)=ustar
                w(3)=pstar
                w(4)=wl(4)
            end if
        end if
    else if (dxdt.lt.lambda(2) .and. dxdt.ge.ustar) then
        if (pstar.gt.wr(3)) then    !right shock
            w(1)=rhostar(2)
            w(2)=ustar
            w(3)=pstar
            w(4)=wr(4)
        else                        !right rarefaction
            csstar=eos_tabulated_sound_speed(rhostar(2),pstar)
            if (dxdt.gt.ustar+csstar) then
                call find_xright_rare_fan(pstar,wr,dxdt,w)
            else
                w(1)=rhostar(2)
                w(2)=ustar
                w(3)=pstar
                w(4)=wr(4)
            end if
        end if
    else
        w=wr
    end if
end subroutine eos_xstate

end module eosriemannsolver
