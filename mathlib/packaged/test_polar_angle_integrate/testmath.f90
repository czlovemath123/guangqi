program testmath
use mathlib
implicit none
real(8), allocatable :: I_bundle(:)
integer :: i,n
real(8) :: mom0,mom1(3)
n=360
allocate(I_bundle(n))
do i=1,n
    if (i.ge.91.and.i.le.271) then
        I_bundle(i)=1d0
    else
        I_bundle(i)=0d0
    end if
end do
call polar_angle_0th_mom_integrate(I_bundle,mom0)
call polar_angle_1st_mom_integrate(I_bundle,mom1)
print *,mom0/pi
print *,mom1/pi

contains

end program testmath
