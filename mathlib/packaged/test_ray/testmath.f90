program testmath
use mathlib
implicit none
real(8) :: s0(3),s1(3),e(3),c(3),ds,mom0,r_forbid,theta,dtheta,dr
real(8), allocatable :: r_depend(:),d(:),I_bundle(:)
integer :: code,n_intersect,n,i,n_polar_rays
logical :: blocked
dr=0.1d0
s0=(/0d0,2d0,0d0/)
c=(/0d0,0d0,0d0/)
r_forbid=1d0
n=20
n_polar_rays=400
allocate(r_depend(n))
do i=1,n
    r_depend=r_forbid+0.08d0*i
end do
dtheta=2d0*pi/n_polar_rays
open(unit=11,file='something.txt',status='replace',action='write')
do i=1,n_polar_rays
    theta=(i-1)*dtheta
    e=(/sin(theta),cos(theta),0d0/)
    call ray_sphere_3d_forbid_depend(s0,e,c,r_forbid,r_depend,s1,dr)
    write(unit=11,fmt='(3ES18.6E2)') s1
end do
close(11)

!do i=1,n
!    r(i)=0.8d0+i*0.1d0
!    call ray_block_single(s,e,c,r(i),ds,blocked)
!    print *, ds,blocked
!end do

!allocate(I_bundle(360))
!I_bundle=1d0
!call polar_angle_0th_mom_integrate(I_bundle,mom0)
!print *, (mom0-4d0*pi)/4d0/pi
!deallocate(I_bundle)

!call ray_sphere_3d_group(s,e,c,r,d,n_intersect)
!print *,n_intersect
!print *,d
deallocate(r_depend)
contains
end program testmath
