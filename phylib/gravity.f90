module gravity
use datastructure
use mathlib
use phylib
use eos
implicit none

contains

subroutine gravity_source(blk)
    !gravity source is second order, t=0 is stored in s_grav1 and t=dt is stored in s_grav2
    type(blockdef), pointer :: blk
    character(len=128) :: alert
    if (igravity==1) then
        !no self gravity
        call gravity_no_self_block(blk)
    else if (igravity==2) then
        !uniform gravity
        call gravity_uniform(blk)
    else
        !not implemented yet
    end if
end subroutine gravity_source

subroutine gravity_no_self_block(blk)
    type(blockdef), pointer :: blk
    if (igeometry==2) then
    else if (igeometry==1) then
    else
    end if
end subroutine gravity_no_self_block

subroutine gravity_uniform(blk)
    type(blockdef), pointer :: blk
    integer :: i,j
    blk%s_grav=0d0
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            blk%s_grav(i,j,1,1)=0d0
            blk%s_grav(i,j,1,2)=0d0
            blk%s_grav(i,j,1,3)=-blk%w(i,j,1,1)*g_uniform
            blk%s_grav(i,j,1,4)=0d0
            blk%s_grav(i,j,1,5)=0d0
        end do
    end do
end subroutine gravity_uniform

subroutine ff_timescale(m,r,dt)
    real(8) :: m,r,dt
    dt=pi/2d0*pow(r,1.5d0)/sqrt(2*gr*m)
end subroutine ff_timescale

end module gravity
