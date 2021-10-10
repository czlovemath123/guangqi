program band
implicit none
integer, parameter :: ku=3,kl=3,nx=512,n=1024,n_example=6,geometry=3
integer, parameter :: transmissive=1,reflective=2,extrapolate=3,specified=9
real(8), parameter :: pi=3.141592652d0,a_rad=7.5646d-15
real(8), dimension(0:nx+1) :: Erad,Frad,feddington,chi,temp
real(8) :: b(n),a(n,n),ab(2*kl+ku+1,n),c_light,alpha,theta(6),phi(6),omega(5),fl,fm,fr,dt,dx,  &
    chim,r,rmin,dr,lscale,vscale,timescale,ler,lfr,r_ratio,rmax,s
integer :: i,j,k,m,ldb,nrhs,ldab,info,ipiv(n),l
integer :: kl_example,ku_example,ipiv_example(n_example),lboundary,rboundary
real(8) :: a_example(6,6),b_example(6),ab_example(7,6)
nrhs=1
ldab=2*kl+ku+1
vscale=1d0
timescale=1d0
lscale=vscale*timescale
rmin=5d12/lscale
rmax=5d13/lscale
ldb=n
c_light=2.998d10/vscale
dt=1d1/timescale
dx=(rmax-rmin)/nx
dr=dx
alpha=dt/dx
!assign physical value
do i=0,nx+1
    feddington(i)=1d0/1d0
    Erad(i)=1d-2
    Frad(i)=Erad(i)   !Erad(i)
    chi(i)=1d-15
    temp(i)=1d4
end do
chi(201)=1d-16
temp(201:210)=1d5
!do i=200,210
!    feddington(i)=1d0-2d0/3d0*(i-200)/10
!end do
!do i=290,300
!    feddington(i)=1d0/3d0+2d0/3d0*(i-290)/10
!end do
!temp(1)=0d0
!temp(nx)=0d0
Erad(0)=1d2
Frad(0)=1d0*Erad(0)/1d0
Erad(nx+1)=1d2
Frad(nx+1)=-Erad(nx+1)
feddington(0)=1d0/1d0
feddington(nx+1)=1d0/1d0
!construct matrix
lboundary=9
rboundary=1
do l=1,1
    do i=1,nx
        s=chi(i)*4*pi*dt*planck_function(temp(i))
        !s=0d0!chi(i)*4*pi*dt*Erad(i)*1d10
        if (i==1) then
            if (lboundary==transmissive) then
            else if (lboundary==reflective) then
            else if (lboundary==specified) then
                theta(1)=theta1(alpha,feddington(0),feddington(1))
                theta(4)=theta4(alpha,feddington(0),feddington(1))
                phi(1)=phi1(alpha,feddington(0),feddington(1))
                phi(4)=phi4(alpha,feddington(0),feddington(1))
                b(1)=Erad(1)-theta(1)*Erad(0)-theta(4)*Frad(0)+s
                b(2)=Frad(1)-phi(1)*Erad(0)-phi(4)*Frad(0)
            end if
        else if (i==nx) then
            if (rboundary==transmissive) then
                b(n-1)=Erad(nx)+s
                b(n)=Frad(nx)
            else if (rboundary==reflective) then
            else if (rboundary==specified) then
                theta(3)=theta3(alpha,feddington(nx),feddington(nx+1))
                theta(6)=theta6(alpha,feddington(nx),feddington(nx+1))
                phi(3)=phi3(alpha,feddington(nx),feddington(nx+1))
                phi(6)=phi6(alpha,feddington(nx),feddington(nx+1))
                b(n-1)=Erad(nx)-theta(3)*Erad(nx+1)-theta(6)*Frad(nx+1)+s
                b(n)=Frad(nx)-phi(3)*Erad(nx+1)-phi(6)*Frad(nx+1)
            end if
        else
            b(2*(i-1)+1)=Erad(i)+s
            b(2*(i-1)+2)=Frad(i)
        end if
    end do
    !write(*,'(20ES12.3E2)') Erad
    !write(*,'(20ES12.3E2)') Frad
    !write(*,'(20ES12.3E2)') feddington

    a=0d0
    do i=1,nx
        r=rmin+dr*i
        if (i==1) then
            if (lboundary==transmissive) then
            else if (lboundary==reflective) then
            else if (lboundary==specified) then
                fl=feddington(0)
                fm=feddington(1)
                fr=feddington(2)
                chim=chi(1)
                theta(2)=theta2(alpha,fl,fm,fr,chim,dt)
                theta(3)=theta3(alpha,fm,fr)
                theta(5)=theta5(alpha,fl,fm,fr,dt,r)
                theta(6)=theta6(alpha,fm,fr)
                phi(2)=phi2(alpha,fl,fm,fr,dt,r)
                phi(3)=phi3(alpha,fm,fr)
                phi(5)=phi5(alpha,fl,fm,fr,chim,dt)
                phi(6)=phi6(alpha,fm,fr)
                a(1,1)=theta(2)
                a(1,2)=theta(5)
                a(1,3)=theta(3)
                a(1,4)=theta(6)
                a(2,1)=phi(2)
                a(2,2)=phi(5)
                a(2,3)=phi(3)
                a(2,4)=phi(6)
            end if
        else if (i==nx) then
            fl=feddington(nx-1)
            fm=feddington(nx)
            fr=feddington(nx+1)
            chim=chi(nx)
            r_ratio=r**2/(r+dr)**2
            if (rboundary==transmissive) then
                theta(1)=theta1(alpha,fl,fm)
                theta(2)=theta2(alpha,fl,fm,fr,chim,dt)
                theta(3)=theta3(alpha,fm,fr)
                theta(4)=theta4(alpha,fl,fm)
                theta(5)=theta5(alpha,fl,fm,fr,dt,r)
                theta(6)=theta6(alpha,fm,fr)
                phi(1)=phi1(alpha,fl,fm)
                phi(2)=phi2(alpha,fl,fm,fr,dt,r)
                phi(3)=phi3(alpha,fm,fr)
                phi(4)=phi4(alpha,fl,fm)
                phi(5)=phi5(alpha,fl,fm,fr,chim,dt)
                phi(6)=phi6(alpha,fm,fr)
                if (geometry==3) then
                    a(n-1,n-3)=theta(1)
                    a(n-1,n-2)=theta(4)
                    a(n-1,n-1)=theta(2)!+r_ratio*theta(3)
                    a(n-1,n)=theta(5)+r_ratio*(theta(6)+theta(3))
                    a(n,n-3)=phi(1)
                    a(n,n-2)=phi(4)
                    a(n,n-1)=phi(2)!+r_ratio*phi(3)
                    a(n,n)=phi(5)+r_ratio*(phi(6)+phi(3))
                else
                    a(n-1,n-3)=theta(1)
                    a(n-1,n-2)=theta(4)
                    a(n-1,n-1)=theta(2)+theta(3)
                    a(n-1,n)=theta(5)+theta(6)!+theta(3)
                    a(n,n-3)=phi(1)
                    a(n,n-2)=phi(4)
                    a(n,n-1)=phi(2)+phi(3)
                    a(n,n)=phi(5)+phi(6)!+phi(3)
                end if
            else if (rboundary==reflective) then
            else if (rboundary==extrapolate) then
                theta(1)=theta1(alpha,fl,fm)
                theta(2)=theta2(alpha,fl,fm,fr,chim,dt)
                theta(3)=theta3(alpha,fm,fr)
                theta(4)=theta4(alpha,fl,fm)
                theta(5)=theta5(alpha,fl,fm,fr,dt,r)
                theta(6)=theta6(alpha,fm,fr)
                phi(1)=phi1(alpha,fl,fm)
                phi(2)=phi2(alpha,fl,fm,fr,dt,r)
                phi(3)=phi3(alpha,fm,fr)
                phi(4)=phi4(alpha,fl,fm)
                phi(5)=phi5(alpha,fl,fm,fr,chim,dt)
                phi(6)=phi6(alpha,fm,fr)
                a(n-1,n-3)=theta(1)
                a(n-1,n-2)=theta(4)
                a(n-1,n-1)=theta(2)+theta(3)+theta(6)
                a(n-1,n)=theta(5)
                a(n,n-3)=phi(1)
                a(n,n-2)=phi(4)
                a(n,n-1)=phi(2)+phi(3)+phi(6)
                a(n,n)=phi(5)
            else if (rboundary==specified) then
                theta(1)=theta1(alpha,fl,fm)
                theta(2)=theta2(alpha,fl,fm,fr,chim,dt)
                theta(4)=theta4(alpha,fl,fm)
                theta(5)=theta5(alpha,fl,fm,fr,dt,r)
                phi(1)=phi1(alpha,fl,fm)
                phi(2)=phi2(alpha,fl,fm,fr,dt,r)
                phi(4)=phi4(alpha,fl,fm)
                phi(5)=phi5(alpha,fl,fm,fr,chim,dt)
                a(n-1,n-3)=theta(1)
                a(n-1,n-2)=theta(4)
                a(n-1,n-1)=theta(2)
                a(n-1,n)=theta(5)
                a(n,n-3)=phi(1)
                a(n,n-2)=phi(4)
                a(n,n-1)=phi(2)
                a(n,n)=phi(5)
            end if
        else
            fl=feddington(i-1)
            fm=feddington(i)
            fr=feddington(i+1)
            chim=chi(i)
            theta(1)=theta1(alpha,fl,fm)
            theta(2)=theta2(alpha,fl,fm,fr,chim,dt)
            theta(3)=theta3(alpha,fm,fr)
            theta(4)=theta4(alpha,fl,fm)
            theta(5)=theta5(alpha,fl,fm,fr,dt,r)
            theta(6)=theta6(alpha,fm,fr)
            phi(1)=phi1(alpha,fl,fm)
            phi(2)=phi2(alpha,fl,fm,fr,dt,r)
            phi(3)=phi3(alpha,fm,fr)
            phi(4)=phi4(alpha,fl,fm)
            phi(5)=phi5(alpha,fl,fm,fr,chim,dt)
            phi(6)=phi6(alpha,fm,fr)
            k=2*(i-1)+1
            a(k,k-2)=theta(1)
            a(k,k-1)=theta(4)
            a(k,k)=theta(2)
            a(k,k+1)=theta(5)
            a(k,k+2)=theta(3)
            a(k,k+3)=theta(6)
            k=2*(i-1)+2
            a(k,k-3)=phi(1)
            a(k,k-2)=phi(4)
            a(k,k-1)=phi(2)
            a(k,k)=phi(5)
            a(k,k+1)=phi(3)
            a(k,k+2)=phi(6)
        end if
    end do
    !print *,'a matrix'
    !write(*,'(20ES10.2E2)') a

    ab=0d0
    do i=1,n
        do j=1,n
            if (i>=max(1,j-ku).and.i<=min(n,j+kl)) then
                ab(kl+ku+1+i-j,j)=a(i,j)
            end if
        end do
    end do
    !print *,'band storage'
    !write(*,'(10ES10.2E2)') ab
    !print *,'rhs'
    !write(*,'(10ES12.3E2)') b
    call dgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
    !print *,info
    do i=1,nx
        Erad(i)=b(2*(i-1)+1)
        Frad(i)=b(2*(i-1)+2)
    end do
    !ler=oscillation(Erad,nx,Erad(0))
    !lfr=oscillation(Frad,nx,Frad(0))
    !print *,ler,lfr
end do
    write(*,'(10ES16.6E2)') Frad    !*c_light
    write(*,'(10ES16.6E2)') Erad



!an example of dgbsv

!b_example=(/1,2,3,4,5,6/)
!do i=1,6
!    do j=1,6
!        if (i==j) then
!            a_example(i,j)=1
!        else if (i==j-1.or.i==j-2) then
!            a_example(i,j)=0.2
!        else if (i==j+1.or.i==j+2) then
!            a_example(i,j)=0.1
!        else
!            a_example(i,j)=0
!        end if
!    end do
!end do
!
!print *,'a matrix'
!write(*,'(6ES18.3E2)') a_example
!
!ab_example=0
!ku_example=2
!kl_example=2
!do i=1,6
!    do j=1,6
!        if (i>=max(1,j-ku_example).and.i<=min(6,j+kl_example)) then
!            ab_example(kl_example+ku_example+1+i-j,j)=a_example(i,j)
!        end if
!    end do
!end do
!
!print *, 'banded storage'
!write(*,'(7ES18.3E2)') ab_example
!
!nrhs=1
!ldb=n_example
!ldab=2*kl_example+ku_example+1
!write(*,'(6ES18.5E3)') b_example
!call dgbsv(n_example,kl_example,ku_example,nrhs,ab_example,ldab,ipiv_example,b_example,ldb,info)
!print *,info
!write(*,'(6ES18.5E3)') b_example

contains

function oscillation(arr,n,v)
    real(8) :: arr(:),v,oscillation,d,l2
    integer :: i,n
    d=1d0/n
    oscillation=0d0
    do i=1,n
        l2=(arr(i)-v)**2*d
        oscillation=oscillation+l2
    end do
end function oscillation

function theta1(alpha,fl,fr)
    real(8) :: alpha,fl,fr,theta1
    theta1=-alpha*c_light*sqrt(fl*fr)/(sqrt(fl)+sqrt(fr))
end function theta1

function theta2(alpha,fl,fm,fr,chi,dt)
    real(8) :: alpha,fl,fm,fr,theta2,chi,dt
    theta2=1d0+alpha*c_light*(sqrt(fl*fm)/(sqrt(fl)+sqrt(fm))+sqrt(fm*fr)/(sqrt(fm)+sqrt(fr)))+chi*c_light*dt
end function theta2

function theta3(alpha,fl,fr)
    real(8) :: alpha,fl,fr,theta3
    theta3=-alpha*c_light*sqrt(fl*fr)/(sqrt(fl)+sqrt(fr))
end function theta3

function theta4(alpha,fl,fr)
    real(8) :: alpha,fl,fr,theta4
    theta4=-alpha*c_light*sqrt(fr)/(sqrt(fl)+sqrt(fr))
end function theta4

function theta5(alpha,fl,fm,fr,dt,r)
    real(8) :: alpha,fl,fm,fr,theta5
    real(8), optional :: dt,r
    theta5=alpha*c_light*(sqrt(fr)/(sqrt(fm)+sqrt(fr))-sqrt(fl)/(sqrt(fl)+sqrt(fm)))
    if (geometry==3) then
        theta5=theta5+2*c_light*dt/r
    end if
end function theta5

function theta6(alpha,fl,fr)
    real(8) :: alpha,fl,fr,theta6
    theta6=alpha*c_light*sqrt(fl)/(sqrt(fl)+sqrt(fr))
end function theta6

function phi1(alpha,fl,fr)
    real(8) :: alpha,fl,fr,phi1
    phi1=-alpha*c_light*fl*sqrt(fr)/(sqrt(fl)+sqrt(fr))
end function phi1

function phi2(alpha,fl,fm,fr,dt,r)
    real(8) :: alpha,fl,fm,fr,phi2
    real(8), optional :: dt,r
    phi2=alpha*c_light*(fm*sqrt(fr)/(sqrt(fm)+sqrt(fr))-fm*sqrt(fl)/(sqrt(fl)+sqrt(fm)))
    if (geometry==3) then
        phi2=phi2+2*fm*c_light*dt/r
    end if
end function phi2

function phi3(alpha,fl,fr)
    real(8) :: alpha,fl,fr,phi3
    phi3=alpha*c_light*fr*sqrt(fl)/(sqrt(fl)+sqrt(fr))
end function phi3

function phi4(alpha,fl,fr)
    real(8) :: alpha,fl,fr,phi4
    phi4=-alpha*c_light*sqrt(fl*fr)/(sqrt(fl)+sqrt(fr))
end function phi4

function phi5(alpha,fl,fm,fr,chi,dt)
    real(8) :: alpha,fl,fm,fr,phi5,chi,dt
    phi5=1d0+alpha*c_light*(sqrt(fl*fm)/(sqrt(fl)+sqrt(fm))+sqrt(fm*fr)/(sqrt(fm)+sqrt(fr)))+chi*c_light*dt
end function phi5

function phi6(alpha,fl,fr)
    real(8) :: alpha,fl,fr,phi6
    phi6=-alpha*c_light*sqrt(fl*fr)/(sqrt(fl)+sqrt(fr))
end function phi6

function planck_function(temp)
    !rad intensity integrated over frequency and assume blackbody
    real(8) :: planck_function,temp
    planck_function=a_rad*c_light/4d0/pi*temp**4
end function planck_function

end program band
