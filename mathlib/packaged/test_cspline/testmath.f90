program testmath
use mathlib
implicit none
real(8), allocatable :: x(:),y(:),t1(:),t2(:),x_cspline(:),x_monotone_cspline(:),  &
    v(:),v_monotone(:),deri(:),deri_monotone(:),v_check(:),deri_check(:),inte(:),  &
    inte_check(:)
real(8) :: dx,v0
integer :: i,j,n,m
n=10
m=9
allocate(x(n),y(n))
do i=1,n
    x(i)=i*1d0-5.5d0
    if (i.le.5) then
        y(i)=1d0
    else
        y(i)=-1d0
    end if
end do
call tangent(x,y,t1)
call monotone_tangent(x,y,t2)
do i=1,n
    print *,t1(i),t2(i)
end do
!cubic spline
call cspline_init(x,m,x_cspline,v,deri,v_check,deri_check)
call cspline_init(x,m,x_cspline,inte,inte_check)
call cspline(x,y,t1,x_cspline,v_check)
do i=1,n-1
    do j=1,n
        call cspline(x(i),x(i+1),y(i),y(i+1),t1(i),t1(i+1),x_cspline((i-1)*n+j),v((i-1)*n+j))
    end do
end do
call cspline(x(n-1),x(n),y(n-1),y(n),t1(n-1),t1(n),x_cspline((n-1)*n+1),v((n-1)*n+1))
!print *,v-v_check
!the first order derivative of cubic spline
call cspline_derivative1(x,y,t1,x_cspline,deri_check)
do i=1,n-1
    do j=1,n
        call cspline_derivative1(x(i),x(i+1),y(i),y(i+1),t1(i),t1(i+1),x_cspline((i-1)*n+j),deri((i-1)*n+j))
    end do
end do
call cspline_derivative1(x(n-1),x(n),y(n-1),y(n),t1(n-1),t1(n),x_cspline((n-1)*n+1),deri((n-1)*n+1))
!print *,deri-deri_check
!monotone cubic spline
call cspline_init(x,m,x_monotone_cspline,v_monotone,deri_monotone)
do i=1,n-1
    do j=1,n
        call cspline(x(i),x(i+1),y(i),y(i+1),t2(i),t2(i+1),  &
            x_monotone_cspline((i-1)*n+j),v_monotone((i-1)*n+j))
    end do
end do
call cspline(x(n-1),x(n),y(n-1),y(n),t2(n-1),t2(n),  &
    x_monotone_cspline((n-1)*n+1),v_monotone((n-1)*n+1))
!the first order derivative of monotone cubic spline
do i=1,n-1
    do j=1,n
        call cspline_derivative1(x(i),x(i+1),y(i),y(i+1),t2(i),t2(i+1),x_cspline((i-1)*n+j),deri_monotone((i-1)*n+j))
    end do
end do
call cspline_derivative1(x(n-1),x(n),y(n-1),y(n),t2(n-1),t2(n),x_cspline((n-1)*n+1),deri_monotone((n-1)*n+1))
!integration
call cspline_integ_array(x,y,t1,x_cspline,inte)
call integrate_profile_directive(x_cspline,v_check,inte_check)
write(*,'(10ES15.6E2)'),inte
write(*,'(10ES15.6E2)'),inte_check
write(*,'(10ES15.6E2)'),inte-inte_check
open(unit=11,file='cspline.txt',status='replace',action='write')
open(unit=12,file='monotone_cspline.txt',status='replace',action='write')
open(unit=13,file='original.txt',status='replace',action='write')
open(unit=14,file='cspline_deri.txt',status='replace',action='write')
open(unit=15,file='cspline_deri_monotone.txt',status='replace',action='write')
do i=1,(n-1)*n+1
    write(unit=11,fmt='(3ES18.6E2)') x_cspline(i),v(i)
    write(unit=12,fmt='(3ES18.6E2)') x_monotone_cspline(i),v_monotone(i)
    write(unit=14,fmt='(3ES18.6E2)') x_cspline(i),deri(i)
    write(unit=15,fmt='(3ES18.6E2)') x_cspline(i),deri_monotone(i)
end do
do i=1,n
    write(unit=13,fmt='(3ES18.6E2)') x(i),y(i)
end do
close(11)
close(12)
close(13)
close(14)
close(15)
deallocate(x,y,x_cspline,x_monotone_cspline,v,v_monotone,deri,deri_monotone)
contains

end program testmath
