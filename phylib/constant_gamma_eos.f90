module constant_gamma_eos
use mathlib
use phylib
implicit none
contains

!convert conserved quantity u to primitive quantity w
subroutine utow(u,w)
    real(8), dimension(5) :: u,w
    real(8) :: rho,vx,vy,vz,egv,p
    !for w, it is rho, vx, vy, vz, p
    rho=u(1)
    vx=u(2)/rho
    vy=u(3)/rho
    vz=u(4)/rho
    egv=u(5)-half*rho*(vx*vx+vy*vy+vz*vz)
    p=egv*(gamma_gas-one)
    w(1)=u(1)
    w(2)=vx
    w(3)=vy
    w(4)=vz
    w(5)=p
end subroutine utow

!convert primitive quantity w to conserved quantity u
subroutine wtou(w,u)
    real(8), dimension(5) :: u,w
    real(8) :: rho,vx,vy,vz,E,p
    !for w, it is rho, vx, vy, vz, p
    rho=w(1)
    vx=w(2)
    vy=w(3)
    vz=w(4)
    p=w(5)
    E=p/(gamma_gas-one)+half*rho*(vx*vx+vy*vy+vz*vz)
    u(1)=rho
    u(2)=rho*vx
    u(3)=rho*vy
    u(4)=rho*vz
    u(5)=E
end subroutine wtou

!convert conserved quantities to source form
subroutine utosource(u,source)
    real(8), dimension(5) :: u,source
    real(8) :: rho,vx,vy,vz
    source(1:4)=u(1:4)
    rho=u(1)
    vx=u(2)/rho
    vy=u(3)/rho
    vz=u(4)/rho
    source(5)=u(5)-half*rho*(vx**2+vy**2+vz**2)
end subroutine utosource

!convert source form to conserved quantities
subroutine sourcetou(source,u)
    real(8), dimension(5) :: u,source
    real(8) :: rho,vx,vy,vz
    u(1:4)=source(1:4)
    rho=source(1)
    vx=source(2)/rho
    vy=source(3)/rho
    vz=source(4)/rho
    u(5)=source(5)+half*rho*(vx**2+vy**2+vz**2)
end subroutine sourcetou

!calculate xflux from w
subroutine wtoflux(w,flux)
    real(8), dimension(5) :: w,flux
    real(8) :: rho,vx,vy,vz,egv,p
    rho=w(1)
    vx=w(2)
    vy=w(3)
    vz=w(4)
    p=w(5)
    egv=p/(gamma_gas-1)
    flux(1)=rho*vx
    flux(2)=rho*vx*vx+p
    flux(3)=rho*vx*vy
    flux(4)=rho*vx*vz
    flux(5)=vx*(0.5d0*rho*(vx*vx+vy*vy+vz*vz)+egv+p)
end subroutine wtoflux

function srhot(rho,t)
    !monatomic gas only
    real(8) :: srhot,rho,t,n
    n=rho/maw/amu
    srhot=n*kb*(1.5d0*log(2d0*pi*maw*amu*kb*t/h_planck**2)-log(n)+2.5d0)/rho
end function srhot

function solvets(s,rho)
    !monatomic gas only, s is the specific entropy
    real(8) :: solvets,s,rho,n
    n=rho/maw/amu
    !print *,s,n,s/n/kb
    !print *,exp(2d0/3d0*(s/n/kb-2.5d0+log(n)))
    !print *,h_planck**2/2d0/pi/maw/amu/kb
    solvets=h_planck**2/2d0/pi/maw/amu/kb*exp(2d0/3d0*(rho*s/n/kb-2.5d0+log(n)))
end function solvets

function idealgas_entropy(rho,temp)
    real(8) :: idealgas_entropy,rho,temp
    !idealgas_entropy=kb/maw/amu*(1.5d0*log(four/three/(gamma_gas-one)*pi*maw*amu*kb*temp/h_planck**2)-log(rho/maw/amu)+2.5d0)
    idealgas_entropy=kb/maw/amu*(1d0/(gamma_gas-1d0)*log(2d0*pi*maw*amu*kb*temp/h_planck**2)-log(rho/maw/amu)+gamma_gas/(gamma_gas-1d0))
end function idealgas_entropy

function idealgas_temp_p(rho,p)
    !calculate temperature from rho, p
    real(8) :: rho,p,idealgas_temp_p
    idealgas_temp_p=maw*amu*p/rho/kb
end function idealgas_temp_p

function idealgas_egmrhop(rho,p)
    !calculate erg/g
    real(8) :: rho,p,idealgas_egmrhop
    idealgas_egmrhop=p/rho/(gamma_gas-one)
end function idealgas_egmrhop

function idealgas_egv(p)
    !calculate erg/cm^3
    real(8) :: p,idealgas_egv
    idealgas_egv=p/(gamma_gas-one)
end function idealgas_egv

function idealgas_temp_egv(rho,egv)
    !calculate temperature from rho, egv
    real(8) :: idealgas_temp_egv,rho,egv
    idealgas_temp_egv=egv*(gamma_gas-one)*maw*amu/rho/kb
end function idealgas_temp_egv

function idealgas_prhot(rho,temp)
    real(8) :: rho,temp,idealgas_prhot
    idealgas_prhot=rho*kb*temp/maw/amu
end function idealgas_prhot

function idealgas_p(egv)
    real(8) :: egv,idealgas_p
    idealgas_p=egv*(gamma_gas-one)
end function idealgas_p

end module constant_gamma_eos
