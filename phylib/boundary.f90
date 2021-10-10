module boundary
use datastructure
use mathlib
use phylib
use communication
use problem
implicit none

contains

subroutine applyboundconds()
    call time_dependent_bound_type()
    call applyboundconds_hydro(boundary_hydro)
    if (iradiation==4) then
        call applyboundconds_rad(boundary_rad)
    end if
end subroutine applyboundconds

subroutine applyboundconds_hydro(sub)
    type(blockdef), pointer :: blk
    procedure(condition_hydro) :: sub
    if (nd==1) then
        call applyboundconds_hydro_1d(sub)
    else if (nd==2) then
        call applyboundconds_hydro_2d(sub)
    end if
end subroutine applyboundconds_hydro

subroutine applyboundconds_hydro_1d(sub)
    type(blockdef), pointer :: blk
    procedure(condition_hydro) :: sub
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        if (associated(blk,blk_head)) then
            call coor1lower_hydro_internal(blk,sub)
        end if
        if (associated(blk,blk_tail)) then
            call coor1upper_hydro_internal(blk,sub)
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end subroutine applyboundconds_hydro_1d

subroutine applyboundconds_hydro_2d(sub)
    !sub is the actual hydro boundary condition
    procedure(condition_hydro) :: sub
    call applyboundconds_hydro_2d_internal(sub)
    if (igeometry==2) then
        call applyboundconds_passive_2d_internal(spherical2d_vphi,boundcond_scalar)
    end if
end subroutine applyboundconds_hydro_2d

!hydro internal communication step

subroutine applyboundconds_hydro_2d_internal(sub)
    type(blockdef), pointer :: blk
    procedure(condition_hydro) :: sub
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        select case (blk%loc_type)
        case (1)
            call xlyu_corner_hydro(blk,sub)
            call applyboundconds_hydro_edge(blk,north,sub)
            call applyboundconds_hydro_edge(blk,west,sub)
        case (2)
            call applyboundconds_hydro_edge(blk,north,sub)
        case (3)
            call xuyu_corner_hydro(blk,sub)
            call applyboundconds_hydro_edge(blk,north,sub)
            call applyboundconds_hydro_edge(blk,east,sub)
        case (4)
            call applyboundconds_hydro_edge(blk,west,sub)
        case (5)
            !inside, do nothing
        case (6)
            call applyboundconds_hydro_edge(blk,east,sub)
        case (7)
            call xlyl_corner_hydro(blk,sub)
            call applyboundconds_hydro_edge(blk,west,sub)
            call applyboundconds_hydro_edge(blk,south,sub)
        case (8)
            call applyboundconds_hydro_edge(blk,south,sub)
        case (9)
            call xuyl_corner_hydro(blk,sub)
            call applyboundconds_hydro_edge(blk,east,sub)
            call applyboundconds_hydro_edge(blk,south,sub)
        case default
        end select
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    nullify(blk)
end subroutine applyboundconds_hydro_2d_internal

subroutine applyboundconds_hydro_edge(blk,direction,sub)
    type(blockdef), pointer :: blk
    procedure(condition_hydro) :: sub
    character(len=16) :: direction
    real(8) :: pos(3),t,w(5),u(5),temp,egv
    integer :: i,j,bound_loc
    character(len=1) :: mirror
    logical :: reflect
    if (direction==east) then
        select case (hydro_bound_type(2))
        case (1)
            do j=1,blk_size_ny
                do i=blk_size_nx+1,blk_size_nx+2
                    blk%w(i,j,1,1:5)=blk%w(blk_size_nx*2+1-i,j,1,1:5)
                end do
            end do
        case (2)
            do j=1,blk_size_ny
                do i=blk_size_nx+1,blk_size_nx+2
                    blk%w(i,j,1,1:5)=blk%w(blk_size_nx*2+1-i,j,1,1:5)
                    blk%w(i,j,1,2)=-blk%w(blk_size_nx*2+1-i,j,1,2)
                end do
            end do
        case (9)
            t=time_sys%t
            bound_loc=2
            do j=-1,blk_size_ny+2
                do i=blk_size_nx+1,blk_size_nx+2
                    pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                    call sub(pos,t,w,u,temp,egv,bound_loc)
                    blk%w(i,j,1,1:5)=w
                end do
            end do
        end select
    else if (direction==south) then
        select case (hydro_bound_type(3))
        case (1)
            do j=-1,0
                do i=1,blk_size_nx
                    blk%w(i,j,1,1:5)=blk%w(i,1-j,1,1:5)
                end do
            end do
        case (2)
            do j=-1,0
                do i=1,blk_size_nx
                    blk%w(i,j,1,1:5)=blk%w(i,1-j,1,1:5)
                    blk%w(i,j,1,3)=-blk%w(i,1-j,1,3)
                end do
            end do
        case (9)
            t=time_sys%t
            bound_loc=3
            do j=-1,0
                do i=-1,blk_size_nx+2
                    pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                    call sub(pos,t,w,u,temp,egv,bound_loc)
                    blk%w(i,j,1,1:5)=w
                end do
            end do
        end select
    else if (direction==west) then
        select case (hydro_bound_type(1))
        case (1)
            do j=1,blk_size_ny
                do i=-1,0
                    blk%w(i,j,1,1:5)=blk%w(1-i,j,1,1:5)
                end do
            end do
        case (2)
            do j=1,blk_size_ny
                do i=-1,0
                    blk%w(i,j,1,1:5)=blk%w(1-i,j,1,1:5)
                    blk%w(i,j,1,2)=-blk%w(1-i,j,1,2)
                end do
            end do
        case (9)
            t=time_sys%t
            bound_loc=1
            do j=-1,blk_size_ny+2
                do i=-1,0
                    pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                    call sub(pos,t,w,u,temp,egv,bound_loc)
                    blk%w(i,j,1,1:5)=w
                end do
            end do
        end select
    else if (direction==north) then
        select case (hydro_bound_type(4))
        case (1)
            do j=blk_size_ny+1,blk_size_ny+2
                do i=1,blk_size_nx
                    blk%w(i,j,1,1:5)=blk%w(i,2*blk_size_ny+1-j,1,1:5)
                end do
            end do
        case (2)
            do j=blk_size_ny+1,blk_size_ny+2
                do i=1,blk_size_nx
                    blk%w(i,j,1,1:5)=blk%w(i,2*blk_size_ny+1-j,1,1:5)
                    blk%w(i,j,1,3)=-blk%w(i,2*blk_size_ny+1-j,1,3)
                end do
            end do
        case (9)
            t=time_sys%t
            bound_loc=4
            do j=blk_size_ny+1,blk_size_ny+2
                do i=-1,blk_size_nx+2
                    pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                    call sub(pos,t,w,u,temp,egv,bound_loc)
                    blk%w(i,j,1,1:5)=w
                end do
            end do
        end select
    end if
end subroutine applyboundconds_hydro_edge

subroutine coor1lower_hydro_internal(blk,sub)
    !the lower boundary condition of the first coordinate
    type(blockdef), pointer :: blk,blk_temp
    procedure(condition_hydro) :: sub
    real(8) :: pos(3),t,w(5),u(5),temp,egv
    integer :: i,j,bound_loc
    character(len=128) :: alert
    character(len=1) :: mirror
    logical :: reflect
    bound_loc=1
    select case(hydro_bound_type(1))
    case(1)     !transmissive on the lower boundary of 1st coordinate
        do i=0,blk_xlb,-1
            blk%w(i,1,1,:)=blk%w(1-i,1,1,:)
            blk%u(i,1,1,:)=blk%u(1-i,1,1,:)
            blk%temp(i,1,1)=blk%temp(1-i,1,1)
            blk%egv(i,1,1)=blk%egv(1-i,1,1)
        end do
    case(2)     !reflective on the lower boundary of 1st coordinate
        do i=0,blk_xlb,-1
            blk%w(i,1,1,:)=blk%w(1-i,1,1,:)
            blk%w(i,1,1,2)=-blk%w(1-i,1,1,2)        !overwrite
            blk%u(i,1,1,:)=blk%u(1-i,1,1,:)
            blk%u(i,1,1,2)=-blk%u(1-i,1,1,2)        !overwrite
            blk%temp(i,1,1)=blk%temp(1-i,1,1)
            blk%egv(i,1,1)=blk%egv(1-i,1,1)
        end do
    case(3)     !periodic, not common
    case(9)     !user specified boundary condition
        t=time_sys%t
        do i=-1,0
            pos=(/blk%x_center(i),0d0,0d0/)
            call sub(pos,t,w,u,temp,egv,bound_loc)
            blk%w(i,1,1,1:5)=w
            blk%u(i,1,1,1:5)=u
            blk%temp(i,1,1)=temp
            blk%egv(i,1,1)=egv
        end do
    end select
end subroutine coor1lower_hydro_internal

subroutine coor1upper_hydro_internal(blk,sub)
    !the upper boundary condition of the first coordinate
    type(blockdef), pointer :: blk,blk_temp
    procedure(condition_hydro) :: sub
    real(8) :: pos(3),t,w(5),u(5),temp,egv
    integer :: ijk(3),i,j,bound_loc
    character(len=128) :: alert
    character(len=1) :: mirror
    logical :: reflect
    bound_loc=2
    select case(hydro_bound_type(2))
    case(1)     !transmissive on the upper boundary of 1st coordinate
        do i=blk_size_nx+1,blk_xub
            blk%w(i,1,1,:)=blk%w(2*blk_size_nx-i+1,1,1,:)
            blk%u(i,1,1,:)=blk%u(2*blk_size_nx-i+1,1,1,:)
            blk%temp(i,1,1)=blk%temp(2*blk_size_nx-i+1,1,1)
            blk%egv(i,1,1)=blk%egv(2*blk_size_nx-i+1,1,1)
        end do
    case(2)     !reflective on the upper boundary of 1st coordinate
        do i=blk_size_nx+1,blk_xub
            blk%w(i,1,1,:)=blk%w(2*blk_size_nx-i+1,1,1,:)
            blk%w(i,1,1,2)=-blk%w(2*blk_size_nx-i+1,1,1,2)          !overwrite
            blk%u(i,1,1,:)=blk%u(2*blk_size_nx-i+1,1,1,:)
            blk%u(i,1,1,2)=-blk%u(2*blk_size_nx-i+1,1,1,2)          !overwrite
            blk%temp(i,1,1)=blk%temp(2*blk_size_nx-i+1,1,1)
            blk%egv(i,1,1)=blk%egv(2*blk_size_nx-i+1,1,1)
        end do
    case(3)     !periodic, not common
    case(9)     !user specified boundary condition
        t=time_sys%t
        do i=blk_size_nx+1,blk_xub
            pos=(/blk%x_center(i),0d0,0d0/)
            call sub(pos,t,w,u,temp,egv,bound_loc)
            blk%w(i,1,1,1:5)=w
            blk%u(i,1,1,1:5)=u
            blk%temp(i,1,1)=temp
            blk%egv(i,1,1)=egv
        end do
    end select
end subroutine coor1upper_hydro_internal

subroutine xlyl_corner_hydro(blk,sub)
    !domain corner, no communication cases
    type(blockdef), pointer :: blk
    procedure(condition_hydro) :: sub
    real(8) :: pos(3),w(5),u(5),temp,egv,t
    real(8), dimension(2,2,5) :: w_in,w_out
    integer :: i,j,bound_loc
    character(len=1) :: mirror
    logical :: reflect
    if (hydro_bound_type(1)==9) then
        t=time_sys%t
        bound_loc=1
        do j=-1,0
            do i=-1,0
                pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                call sub(pos,t,w,u,temp,egv,bound_loc)
                blk%w(i,j,1,1:5)=w
            end do
        end do
    else if (hydro_bound_type(3)==9) then
        t=time_sys%t
        bound_loc=3
        do j=-1,0
            do i=-1,0
                pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                call sub(pos,t,w,u,temp,egv,bound_loc)
                blk%w(i,j,1,1:5)=w
            end do
        end do
    else
        call extract_corner_hydro(blk,sw,w_in)
        if (hydro_bound_type(1)==1.and.hydro_bound_type(3)==1) then
            mirror='x';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(1)==1.and.hydro_bound_type(3)==2) then
            mirror='x';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(1)==2.and.hydro_bound_type(3)==1) then
            mirror='x';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(1)==2.and.hydro_bound_type(3)==2) then
            mirror='x';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        end if
    end if
end subroutine xlyl_corner_hydro

subroutine xuyl_corner_hydro(blk,sub)
    !domain corner, no communication cases
    type(blockdef), pointer :: blk
    procedure(condition_hydro) :: sub
    real(8) :: pos(3),w(5),u(5),temp,egv,t
    real(8), dimension(2,2,5) :: w_in,w_out
    integer :: i,j,bound_loc
    character(len=1) :: mirror
    logical :: reflect
    if (hydro_bound_type(2)==9) then
        t=time_sys%t
        bound_loc=2
        do j=-1,0
            do i=blk_size_nx+1,blk_size_nx+2
                pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                call sub(pos,t,w,u,temp,egv,bound_loc)
                blk%w(i,j,1,1:5)=w
            end do
        end do
    else if (hydro_bound_type(3)==9) then
        t=time_sys%t
        bound_loc=3
        do j=-1,0
            do i=blk_size_nx+1,blk_size_nx+2
                pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                call sub(pos,t,w,u,temp,egv,bound_loc)
                blk%w(i,j,1,1:5)=w
            end do
        end do
    else
        call extract_corner_hydro(blk,se,w_in)
        if (hydro_bound_type(2)==1.and.hydro_bound_type(3)==1) then
            mirror='x';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(2)==1.and.hydro_bound_type(3)==2) then
            mirror='x';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(2)==2.and.hydro_bound_type(3)==1) then
            mirror='x';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(2)==2.and.hydro_bound_type(3)==2) then
            mirror='x';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        end if
    end if
end subroutine xuyl_corner_hydro

subroutine xlyu_corner_hydro(blk,sub)
    !domain corner, no communication cases
    type(blockdef), pointer :: blk
    procedure(condition_hydro) :: sub
    real(8) :: pos(3),w(5),u(5),temp,egv,t
    real(8), dimension(2,2,5) :: w_in,w_out
    integer :: i,j,bound_loc
    character(len=1) :: mirror
    logical :: reflect
    if (hydro_bound_type(1)==9) then
        t=time_sys%t
        bound_loc=1
        do j=blk_size_ny+1,blk_size_ny+2
            do i=-1,0
                pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                call sub(pos,t,w,u,temp,egv,bound_loc)
                blk%w(i,j,1,1:5)=w
            end do
        end do
    else if (hydro_bound_type(4)==9) then
        t=time_sys%t
        bound_loc=4
        do j=blk_size_ny+1,blk_size_ny+2
            do i=-1,0
                pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                call sub(pos,t,w,u,temp,egv,bound_loc)
                blk%w(i,j,1,1:5)=w
            end do
        end do
    else
        call extract_corner_hydro(blk,nw,w_in)
        if (hydro_bound_type(1)==1.and.hydro_bound_type(4)==1) then
            mirror='x';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(1)==1.and.hydro_bound_type(4)==2) then
            mirror='x';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(1)==2.and.hydro_bound_type(4)==1) then
            mirror='x';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(1)==2.and.hydro_bound_type(4)==2) then
            mirror='x';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        end if
    end if
end subroutine xlyu_corner_hydro

subroutine xuyu_corner_hydro(blk,sub)
    !domain corner, no communication cases
    type(blockdef), pointer :: blk
    procedure(condition_hydro) :: sub
    real(8) :: pos(3),w(5),u(5),temp,egv,t
    real(8), dimension(2,2,5) :: w_in,w_out
    integer :: i,j,bound_loc
    character(len=1) :: mirror
    logical :: reflect
    if (hydro_bound_type(2)==9) then
        t=time_sys%t
        bound_loc=2
        do j=blk_size_ny+1,blk_size_ny+2
            do i=blk_size_nx+1,blk_size_nx+2
                pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                call sub(pos,t,w,u,temp,egv,bound_loc)
                blk%w(i,j,1,1:5)=w
            end do
        end do
    else if (hydro_bound_type(4)==9) then
        t=time_sys%t
        bound_loc=4
        do j=blk_size_ny+1,blk_size_ny+2
            do i=blk_size_nx+1,blk_size_nx+2
                pos=(/blk%x_center(i),blk%y_center(j),0d0/)
                call sub(pos,t,w,u,temp,egv,bound_loc)
                blk%w(i,j,1,1:5)=w
            end do
        end do
    else
        call extract_corner_hydro(blk,ne,w_in)
        if (hydro_bound_type(2)==1.and.hydro_bound_type(4)==1) then
            mirror='x';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(2)==1.and.hydro_bound_type(4)==2) then
            mirror='x';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(2)==2.and.hydro_bound_type(4)==1) then
            mirror='x';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.false.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        else if (hydro_bound_type(2)==2.and.hydro_bound_type(4)==2) then
            mirror='x';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
            w_in=w_out
            mirror='y';reflect=.true.
            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
        end if
    end if
end subroutine xuyu_corner_hydro

subroutine applyboundconds_rad(sub)
    procedure(condition_rad) :: sub
    if (nd==1) then
        call applyboundconds_rad_1d(sub)
    else if (nd==2) then
        call applyboundconds_rad_2d(sub)
    end if
end subroutine applyboundconds_rad

subroutine applyboundconds_rad_1d(sub)
    type(blockdef), pointer :: blk
    procedure(condition_rad) :: sub
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        if (associated(blk,blk_head)) then
            call coor1lower_rad_boundary(blk,sub)
        end if
        if (associated(blk,blk_tail)) then
            call coor1upper_rad_boundary(blk,sub)
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end subroutine applyboundconds_rad_1d

subroutine applyboundconds_rad_2d(sub)
    type(blockdef), pointer :: blk
    procedure(condition_rad) :: sub
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        select case (blk%loc_type)
        case (1)
        case (2)
        case (3)
        case (4)
        case (5)
        case (6)
        case (7)
        case (8)
        case (9)
        case default
        end select
    end do
end subroutine applyboundconds_rad_2d

!radiation transfer usually use implicit solver, their boundary conditions are implicitly applied in their solver
!the boundary conditions below only apply when there is an explicit boundary condition
subroutine coor1lower_rad_boundary(blk,sub)
    type(blockdef), pointer :: blk
    procedure(condition_rad) :: sub
    real(8) :: pos(3),t,Frad(3),feddington
    real(8), dimension(:), allocatable :: Erad
    integer :: ijk(3),i,j
    character(len=128) :: alert
    if (nd==1) then
        select case(rad_bound_type(1))
        case(9)     !user specified boundary condition
            t=time_sys%t
            pos=(/blk%x_center(0),0d0,0d0/)
            if (iradiation==1) then
                allocate(Erad(1))
                call sub(pos,t,Erad,Frad,feddington)
                blk%Erad(0,1,1)=Erad(1)
                deallocate(Erad)
            else if (iradiation==3) then
            else if (iradiation==4) then
                allocate(Erad(1))
                call sub(pos,t,Erad,Frad,feddington)
                blk%Erad(0,1,1)=Erad(1)
                deallocate(Erad)
            else if (iradiation==5) then
            end if
        case(1)     !transmissive boundary condition
            alert="1d transimissive lower radiation boundary not implemented"
            call abort_achilles(alert)
        case(2)     !reflective boundary condition
            if (iradiation==1.or.iradiation==4) then
                blk%Erad(0,1,1)=blk%Erad(1,1,1)
            else
            end if
        case default
        end select
    else if (nd==2) then
        select case(rad_bound_type(1))
        case(9)     !user specifiec boundary condition
            t=time_sys%t
            do j=0,ny+1
                pos=(/blk%x_center(0),blk%y_center(j),0d0/)
                call sub(pos,t,Erad,Frad,feddington)
            end do
        case default
            !the boundary condition is implicit
        end select
    else
    end if
end subroutine coor1lower_rad_boundary

subroutine coor1upper_rad_boundary(blk,sub)
    type(blockdef), pointer :: blk
    procedure(condition_rad) :: sub
    real(8) :: pos(3),t,Frad(3),feddington
    real(8), dimension(:), allocatable :: Erad
    integer :: ijk(3),i,j
    character(len=128) :: alert
    if (nd==1) then
        select case(rad_bound_type(2))
        case(9)     !user specified boundary condition
            t=time_sys%t
            pos=(/blk%x_center(blk_size_nx+1),0d0,0d0/)
            if (iradiation==1) then
                allocate(Erad(1))
                call sub(pos,t,Erad,Frad,feddington)
                blk%Erad(blk_size_nx+1,1,1)=Erad(1)
                deallocate(Erad)
            else if (iradiation==3) then
            else if (iradiation==4) then
                allocate(Erad(1))
                call sub(pos,t,Erad,Frad,feddington)
                blk%Erad(blk_size_nx+1,1,1)=Erad(1)
                deallocate(Erad)
            else if (iradiation==5) then
            end if
        case(1)     !transmissive boundary condition
            !defined implicitly in radiation_common_functions.f90 calculate_fld_conductivity
        case(2)     !reflective boundary condition
            if (iradiation==5) then
                !blk%Erad_mg(:,blk_size_nx+1,1,1)=blk%Erad_mg(:,blk_size_nx,1,1)
            else
                blk%Erad(blk_size_nx+1,1,1)=blk%Erad(blk_size_nx,1,1)
            end if
        case default
        end select
    else if (nd==2) then
        select case(rad_bound_type(1))
        case(9)     !user specifiec boundary condition
            t=time_sys%t
            do j=0,ny+1
                pos=(/blk%x_center(blk_size_nx+1),blk%y_center(j),0d0/)
                call sub(pos,t,Erad,Frad,feddington)
            end do
        case default
            !the boundary condition is implicit
        end select
    else
    end if
end subroutine coor1upper_rad_boundary

subroutine coor2lower_rad_boundary(blk,sub)
    !only explicit boundary condition
    type(blockdef), pointer :: blk
    procedure(condition_rad) :: sub
    real(8) :: pos(3),t,Frad(3),feddington
    real(8), dimension(:), allocatable :: Erad
    integer :: ijk(3),i,j
    character(len=128) :: alert
    if (nd==1) then
        alert='no 2nd coordinate for 1d sim'
        call abort_achilles(alert)
    else if (nd==2) then
        select case(rad_bound_type(3))
        case(9)     !user specifiec boundary condition
            t=time_sys%t
            do i=0,nx+1
                pos=(/blk%x_center(i),blk%y_center(0),0d0/)
                if (iradiation==1) then
                    allocate(Erad(1))
                    call sub(pos,t,Erad,Frad,feddington)
                    blk%Erad(i,0,1)=Erad(1)
                    deallocate(Erad)
                else if (iradiation==3) then
                end if
            end do
        case default
            !the boundary condition is implicit
        end select
    else
    end if
end subroutine coor2lower_rad_boundary

subroutine coor2upper_rad_boundary(blk,sub)
    !only explicit boundary condition
    type(blockdef), pointer :: blk
    procedure(condition_rad) :: sub
    real(8) :: pos(3),t,Frad(3),feddington
    real(8), dimension(:), allocatable :: Erad
    integer :: ijk(3),i,j
    character(len=128) :: alert
    if (nd==1) then
        alert='no 2nd coordinate for 1d sim'
        call abort_achilles(alert)
    else if (nd==2) then
        select case(rad_bound_type(4))
        case(9)     !user specifiec boundary condition
            t=time_sys%t
            do i=0,nx+1
                pos=(/blk%x_center(i),blk%y_center(blk_size_ny+1),0d0/)
                if (iradiation==1) then
                    allocate(Erad(1))
                    call sub(pos,t,Erad,Frad,feddington)
                    blk%Erad(i,ny+1,1)=Erad(1)
                    deallocate(Erad)
                else if (iradiation==3) then
                end if
            end do
        case default
            !the boundary condition is implicit
        end select
    else
    end if
end subroutine coor2upper_rad_boundary

subroutine applyboundconds_passive_2d_internal(sub,sub1)
    !sub1 include all the specified boundary conditions.
    type(blockdef), pointer :: blk
    procedure(condition_scalar), pointer :: sub1
    procedure(extractarray) :: sub
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        select case (blk%loc_type)
        case (1)
            call coor1lower_passive_internal(blk,sub,sub1)
            call coor2upper_passive_internal(blk,sub,sub1)
            call xlyu_corner_passive_internal(blk,sub,sub1)
        case (2)
            call coor2upper_passive_internal(blk,sub,sub1)
        case (3)
            call coor1upper_passive_internal(blk,sub,sub1)
            call coor2upper_passive_internal(blk,sub,sub1)
            call xuyu_corner_passive_internal(blk,sub,sub1)
        case (4)
            call coor1lower_passive_internal(blk,sub,sub1)
        case (5)
            !do nothing
        case (6)
            call coor1upper_passive_internal(blk,sub,sub1)
        case (7)
            call coor1lower_passive_internal(blk,sub,sub1)
            call coor2lower_passive_internal(blk,sub,sub1)
            call xlyl_corner_passive_internal(blk,sub,sub1)
        case (8)
            call coor2lower_passive_internal(blk,sub,sub1)
        case (9)
            call coor1upper_passive_internal(blk,sub,sub1)
            call coor2lower_passive_internal(blk,sub,sub1)
            call xuyl_corner_passive_internal(blk,sub,sub1)
        case default
        end select
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end subroutine applyboundconds_passive_2d_internal

subroutine coor1lower_passive_internal(blk,sub,sub1)
    type(blockdef), pointer :: blk,blk_l,blk_u
    real(8), dimension(:,:,:), pointer :: ap,ap_l,ap_u
    procedure(extractarray) :: sub
    procedure(condition_scalar), pointer :: sub1
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t,scalar
    if (nd==1) then
    else if (nd==2) then
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        if (hydro_bound_type(1)==9) then
            blk_dxyz=blk%dxyz
            blk_coords=blk%pos
            t=time_sys%t
            call sub(blk,ap)
            do j=blk_ylb,blk_yub
                do i=-1,0
                    ijk=(/i,j,1/)
                    call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                    call sub1(pos,t,scalar)
                    ap(i,j,1)=scalar
                end do
            end do
        else if (hydro_bound_type(1)==1.or.hydro_bound_type(1)==2) then
            !allocate(ap_l(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
            !allocate(ap_u(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
            do j=1,blk_size_ny
                do i=-1,0
                    ap(i,j,1)=ap(1-i,j,1)
                end do
            end do
            if (blk%loc_type==1) then
                if (blk%nb_sw%cpu_rank==rank) then
                    blk_l=>blk%nb_sw%blk
                    call sub(blk_l,ap_l)
                end if
            else if (blk%loc_type==7) then
                if (blk%nb_nw%cpu_rank==rank) then
                    blk_u=>blk%nb_nw%blk
                    call sub(blk_u,ap_u)
                end if
            else if (blk%loc_type==4) then
                if (blk%nb_sw%cpu_rank==rank) then
                    blk_l=>blk%nb_sw%blk
                    call sub(blk_l,ap_l)
                end if
                if (blk%nb_nw%cpu_rank==rank) then
                    blk_u=>blk%nb_nw%blk
                    call sub(blk_u,ap_u)
                end if
            end if
            if (associated(blk_l)) then
                do j=-1,0
                    do i=-1,0
                        ap(i,j,1)=ap_l(1-i,blk_size_ny+j,1)
                    end do
                end do
            end if
            if (associated(blk_u)) then
                do j=blk_size_ny+1,blk_size_ny+2
                    do i=-1,0
                        ap(i,j,1)=ap_u(1-i,j-blk_size_ny,1)
                    end do
                end do
            end if
        else if (hydro_bound_type(1)==3) then
            !not common
        end if
        nullify(ap,ap_l,ap_u,blk_l,blk_u)
    end if
end subroutine coor1lower_passive_internal

subroutine coor1upper_passive_internal(blk,sub,sub1)
    type(blockdef), pointer :: blk,blk_l,blk_u
    real(8), dimension(:,:,:), pointer :: ap,ap_l,ap_u
    procedure(extractarray) :: sub
    procedure(condition_scalar), pointer :: sub1
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t,scalar
    if (nd==1) then
    else if (nd==2) then
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        if (hydro_bound_type(2)==9) then
            blk_dxyz=blk%dxyz
            blk_coords=blk%pos
            t=time_sys%t
            call sub(blk,ap)
            do j=blk_ylb,blk_yub
                do i=blk_size_nx+1,blk_size_nx+2
                    ijk=(/i,j,1/)
                    call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                    call sub1(pos,t,scalar)
                    ap(i,j,1)=scalar
                end do
            end do
        else if (hydro_bound_type(2)==1.or.hydro_bound_type(2)==2) then
        end if
        nullify(ap)
    end if
end subroutine coor1upper_passive_internal

subroutine coor2lower_passive_internal(blk,sub,sub1)
    type(blockdef), pointer :: blk,blk_l,blk_u
    real(8), dimension(:,:,:), pointer :: ap,ap_l,ap_u
    procedure(extractarray) :: sub
    procedure(condition_scalar), pointer :: sub1
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t,scalar
    if (nd==1) then
    else if (nd==2) then
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        if (hydro_bound_type(3)==9) then
            blk_dxyz=blk%dxyz
            blk_coords=blk%pos
            t=time_sys%t
            call sub(blk,ap)
            do j=-1,0
                do i=blk_xlb,blk_xub
                    ijk=(/i,j,1/)
                    call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                    call sub1(pos,t,scalar)
                    ap(i,j,1)=scalar
                end do
            end do
        else if (hydro_bound_type(3)==1.or.hydro_bound_type(3)==2) then
        end if
        nullify(ap)
    end if
end subroutine coor2lower_passive_internal

subroutine coor2upper_passive_internal(blk,sub,sub1)
    type(blockdef), pointer :: blk,blk_l,blk_u
    real(8), dimension(:,:,:), pointer :: ap,ap_l,ap_u
    procedure(extractarray) :: sub
    procedure(condition_scalar), pointer :: sub1
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t,scalar
    if (nd==1) then
    else if (nd==2) then
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        if (hydro_bound_type(4)==9) then
            blk_dxyz=blk%dxyz
            blk_coords=blk%pos
            t=time_sys%t
            call sub(blk,ap)
            do j=blk_size_ny+1,blk_yub
                do i=blk_xlb,blk_xub
                    ijk=(/i,j,1/)
                    call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                    call sub1(pos,t,scalar)
                    ap(i,j,1)=scalar
                end do
            end do
        else if (hydro_bound_type(4)==1.or.hydro_bound_type(4)==2) then
        end if
        nullify(ap)
    end if
end subroutine coor2upper_passive_internal

subroutine xlyl_corner_passive_internal(blk,sub,sub1)
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
    procedure(condition_scalar), pointer :: sub1
    procedure(extractarray) :: sub
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t,scalar
    !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    if (hydro_bound_type(1)==9.or.hydro_bound_type(3)==9) then
        blk_dxyz=blk%dxyz
        blk_coords=blk%pos
        t=time_sys%t
        call sub(blk,ap)
        do j=-1,0
            do i=-1,0
                ijk=(/i,j,1/)
                call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                call sub1(pos,t,scalar)
                ap(i,j,1)=scalar
            end do
        end do
    else if (hydro_bound_type(1)==3.and.hydro_bound_type(3)==3) then
        !both periodic boundary, not common
    else if (hydro_bound_type(1)==3.and.hydro_bound_type(3)/=3) then
    else if (hydro_bound_type(1)/=3.and.hydro_bound_type(3)==3) then
        !cylindrical coordinate use this boundary
    else
        call sub(blk,ap)
        do j=-1,0
            do i=-1,0
                ap(i,j,1)=ap(1-i,1-j,1)
            end do
        end do
    end if
    nullify(ap)
end subroutine xlyl_corner_passive_internal

subroutine xuyl_corner_passive_internal(blk,sub,sub1)
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
    procedure(condition_scalar), pointer :: sub1
    procedure(extractarray) :: sub
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t,scalar
    !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    if (hydro_bound_type(2)==9.or.hydro_bound_type(3)==9) then
        blk_dxyz=blk%dxyz
        blk_coords=blk%pos
        t=time_sys%t
        call sub(blk,ap)
        do j=-1,0
            do i=blk_size_nx+1,blk_size_nx+2
                ijk=(/i,j,1/)
                call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                call sub1(pos,t,scalar)
                ap(i,j,1)=scalar
            end do
        end do
    else if (hydro_bound_type(2)==3.and.hydro_bound_type(3)==3) then
        !both periodic boundary, not common
    else if (hydro_bound_type(2)==3.and.hydro_bound_type(3)/=3) then
    else if (hydro_bound_type(2)/=3.and.hydro_bound_type(3)==3) then
        !cylindrical coordinate use this boundary
    else
        call sub(blk,ap)
        do j=-1,0
            do i=blk_size_nx+1,blk_size_nx+2
                ap(i,j,1)=ap(2*blk_size_nx-i+1,1-j,1)
            end do
        end do
    end if
    nullify(ap)
end subroutine xuyl_corner_passive_internal

subroutine xlyu_corner_passive_internal(blk,sub,sub1)
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
    procedure(condition_scalar), pointer :: sub1
    procedure(extractarray) :: sub
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t,scalar
    !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    if (hydro_bound_type(1)==9.and.hydro_bound_type(4)==9) then
        blk_dxyz=blk%dxyz
        blk_coords=blk%pos
        t=time_sys%t
        call sub(blk,ap)
        do j=blk_size_ny+1,blk_size_ny+2
            do i=-1,0
                ijk=(/i,j,1/)
                call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                call sub1(pos,t,scalar)
                ap(i,j,1)=scalar
            end do
        end do
    else if (hydro_bound_type(1)==3.and.hydro_bound_type(4)==3) then
        !both periodic boundary, not common
    else if (hydro_bound_type(1)==3.and.hydro_bound_type(4)/=3) then
    else if (hydro_bound_type(1)/=3.and.hydro_bound_type(4)==3) then
        !cylindrical coordinate use this boundary
    else
    end if
    nullify(ap)
end subroutine xlyu_corner_passive_internal

subroutine xuyu_corner_passive_internal(blk,sub,sub1)
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
    procedure(condition_scalar), pointer :: sub1
    procedure(extractarray) :: sub
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t,scalar
    !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    if (hydro_bound_type(2)==9.and.hydro_bound_type(4)==9) then
        blk_dxyz=blk%dxyz
        blk_coords=blk%pos
        t=time_sys%t
        call sub(blk,ap)
        do j=blk_size_ny+1,blk_size_ny+2
            do i=blk_size_nx+1,blk_size_nx+2
                ijk=(/i,j,1/)
                call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                call sub1(pos,t,scalar)
                ap(i,j,1)=scalar
            end do
        end do
    else if (hydro_bound_type(2)==3.and.hydro_bound_type(4)==3) then
        !both periodic boundary, not common
    else if (hydro_bound_type(2)==3.and.hydro_bound_type(4)/=3) then
    else if (hydro_bound_type(2)/=3.and.hydro_bound_type(4)==3) then
        !cylindrical coordinate use this boundary
    else
    end if
    nullify(ap)
end subroutine xuyu_corner_passive_internal

subroutine applyboundconds_passive_2d_edge_external(sub)
    !specify the passive boundary condition, including possible communications
    type(blockdef), pointer :: blk
    type(blockneighbour), pointer :: blk_nb
    procedure(extractarray) :: sub
    integer :: i,loc_type,neighbour_rank,ierr,tag
    real(8), dimension(:,:,:), pointer :: ap
    blk=>blk_processor_head
    if (hydro_bound_type(1)==3.and.hydro_bound_type(2)==3) then
        !not common
    else if (hydro_bound_type(3)==3.and.hydro_bound_type(4)==3) then
        !cylindrical coordinate uses this boundary condition
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        do i=1,np_nblk(rank+1)
            call sub(blk,ap)
            select case (blk%loc_type)
            case (1)
                call communicate_passive_edge_external_send(blk,ap,north)
                call communicate_passive_corner_external_send(blk,ap,ne)
                call coor1lower_passive_external_send(blk,sub)
            case (2)
                call communicate_passive_edge_external_send(blk,ap,north)
                call communicate_passive_corner_external_send(blk,ap,ne)
                call communicate_passive_corner_external_send(blk,ap,nw)
            case (3)
                call communicate_passive_edge_external_send(blk,ap,north)
                call communicate_passive_corner_external_send(blk,ap,nw)
                call coor1upper_passive_external_send(blk,sub)
            case (4)
                call coor1lower_passive_external_send(blk,sub)
            case (5)
                !inside do nothing
            case (6)
                call coor1upper_passive_external_send(blk,sub)
            case (7)
                call communicate_passive_edge_external_send(blk,ap,south)
                call communicate_passive_corner_external_send(blk,ap,se)
                call coor1lower_passive_external_send(blk,sub)
            case (8)
                call communicate_passive_edge_external_send(blk,ap,south)
                call communicate_passive_corner_external_send(blk,ap,se)
                call communicate_passive_corner_external_send(blk,ap,sw)
            case (9)
                call communicate_passive_edge_external_send(blk,ap,south)
                call communicate_passive_corner_external_send(blk,ap,sw)
                call coor1upper_passive_external_send(blk,sub)
            case default
            end select
            if (i/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            call sub(blk,ap)
            select case (blk%loc_type)
            case (1)
                call communicate_passive_edge_external_recv(blk,ap,north)
                call communicate_passive_corner_external_recv(blk,ap,ne)
                call coor1lower_passive_external_recv(blk,sub)
            case (2)
                call communicate_passive_edge_external_recv(blk,ap,north)
                call communicate_passive_corner_external_recv(blk,ap,ne)
                call communicate_passive_corner_external_recv(blk,ap,nw)
            case (3)
                call communicate_passive_edge_external_recv(blk,ap,north)
                call communicate_passive_corner_external_recv(blk,ap,nw)
                call coor1upper_passive_external_recv(blk,sub)
            case (4)
                call coor1lower_passive_external_recv(blk,sub)
            case (5)
                !inside do nothing
            case (6)
                call coor1upper_passive_external_recv(blk,sub)
            case (7)
                call communicate_passive_edge_external_recv(blk,ap,south)
                call communicate_passive_corner_external_recv(blk,ap,se)
                call coor1lower_passive_external_recv(blk,sub)
            case (8)
                call communicate_passive_edge_external_recv(blk,ap,south)
                call communicate_passive_corner_external_recv(blk,ap,se)
                call communicate_passive_corner_external_recv(blk,ap,sw)
            case (9)
                call communicate_passive_edge_external_recv(blk,ap,south)
                call communicate_passive_corner_external_recv(blk,ap,sw)
                call coor1upper_passive_external_recv(blk,sub)
            case default
            end select
            if (i/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
        nullify(ap)
    else
        do i=1,np_nblk(rank+1)
            select case (blk%loc_type)
            case (1)
                call coor1lower_passive_external_send(blk,sub)
                call coor2upper_passive_external_send(blk,sub)
            case (2)
                call coor2upper_passive_external_send(blk,sub)
            case (3)
                call coor1upper_passive_external_send(blk,sub)
                call coor2upper_passive_external_send(blk,sub)
            case (4)
                call coor1lower_passive_external_send(blk,sub)
            case (5)
                !inside, do nothing
            case (6)
                call coor1upper_passive_external_send(blk,sub)
            case (7)
                call coor1lower_passive_external_send(blk,sub)
                call coor2lower_passive_external_send(blk,sub)
            case (8)
                call coor2lower_passive_external_send(blk,sub)
            case (9)
                call coor1upper_passive_external_send(blk,sub)
                call coor2lower_passive_external_send(blk,sub)
            case default
            end select
            if (i/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            select case (blk%loc_type)
            case (1)
                call coor1lower_passive_external_recv(blk,sub)
                call coor2upper_passive_external_recv(blk,sub)
            case (2)
                call coor2upper_passive_external_recv(blk,sub)
            case (3)
                call coor1upper_passive_external_recv(blk,sub)
                call coor2upper_passive_external_recv(blk,sub)
            case (4)
                call coor1lower_passive_external_recv(blk,sub)
            case (5)
                !inside, do nothing
            case (6)
                call coor1upper_passive_external_recv(blk,sub)
            case (7)
                call coor1lower_passive_external_recv(blk,sub)
                call coor2lower_passive_external_recv(blk,sub)
            case (8)
                call coor2lower_passive_external_recv(blk,sub)
            case (9)
                call coor1upper_passive_external_recv(blk,sub)
                call coor2lower_passive_external_recv(blk,sub)
            case default
            end select
            if (i/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
    end if
    nullify(blk)
end subroutine applyboundconds_passive_2d_edge_external

subroutine coor1lower_passive_external_send(blk,sub)
    type(blockdef), pointer :: blk,blk_l,blk_u
    real(8), dimension(:,:,:), pointer :: ap,ap_l,ap_u
    procedure(extractarray) :: sub
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t
    if (nd==1) then
    else if (nd==2) then
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        if (hydro_bound_type(1)==9) then
            !specified in the internal step
        else if (hydro_bound_type(1)==1.or.hydro_bound_type(1)==2) then
            if (blk%loc_type==1) then
            else if (blk%loc_type==7) then
            else if (blk%loc_type==4) then
            end if
        else if (hydro_bound_type(1)==3) then
            !not common
        end if
        nullify(ap)
    end if
end subroutine coor1lower_passive_external_send

subroutine coor1upper_passive_external_send(blk,sub)
    type(blockdef), pointer :: blk,blk_l,blk_u
    real(8), dimension(:,:,:), pointer :: ap,ap_l,ap_u
    procedure(extractarray) :: sub
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t
    if (nd==1) then
    else if (nd==2) then
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        if (hydro_bound_type(2)==9) then
            !specified in the internal step
        else if (hydro_bound_type(2)==1.or.hydro_bound_type(2)==2) then
            if (blk%loc_type==1) then
            else if (blk%loc_type==7) then
            else if (blk%loc_type==4) then
            end if
        else if (hydro_bound_type(2)==3) then
            !not common
        end if
        nullify(ap)
    end if
end subroutine coor1upper_passive_external_send

subroutine coor2lower_passive_external_send(blk,sub)
    type(blockdef), pointer :: blk,blk_l,blk_u
    real(8), dimension(:,:,:), pointer :: ap,ap_l,ap_u
    procedure(extractarray) :: sub
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t
    if (nd==1) then
    else if (nd==2) then
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        if (hydro_bound_type(3)==9) then
            !specified in the internal step
        else if (hydro_bound_type(3)==1.or.hydro_bound_type(3)==2) then
            if (blk%loc_type==1) then
            else if (blk%loc_type==7) then
            else if (blk%loc_type==4) then
            end if
        else if (hydro_bound_type(3)==3) then
            !cylindrical coordinate uses this boundary condition
        end if
        nullify(ap)
    end if
end subroutine coor2lower_passive_external_send

subroutine coor2upper_passive_external_send(blk,sub)
    type(blockdef), pointer :: blk,blk_l,blk_u
    real(8), dimension(:,:,:), pointer :: ap,ap_l,ap_u
    procedure(extractarray) :: sub
    integer :: i,j,ijk(3)
    real(8) :: blk_dxyz(3),blk_coords(3),pos(3),t
    if (nd==1) then
    else if (nd==2) then
        !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
        if (hydro_bound_type(4)==9) then
            !specified in the internal step
        else if (hydro_bound_type(4)==1.or.hydro_bound_type(4)==2) then
            if (blk%loc_type==1) then
            else if (blk%loc_type==7) then
            else if (blk%loc_type==4) then
            end if
        else if (hydro_bound_type(4)==3) then
            !cylindrical coordinate uses this boundary condition
        end if
        nullify(ap)
    end if
end subroutine coor2upper_passive_external_send

subroutine coor1lower_passive_external_recv(blk,sub)
    type(blockdef), pointer :: blk
    procedure(extractarray) :: sub
end subroutine coor1lower_passive_external_recv

subroutine coor1upper_passive_external_recv(blk,sub)
    type(blockdef), pointer :: blk
    procedure(extractarray) :: sub
end subroutine coor1upper_passive_external_recv

subroutine coor2lower_passive_external_recv(blk,sub)
    type(blockdef), pointer :: blk
    procedure(extractarray) :: sub
end subroutine coor2lower_passive_external_recv

subroutine coor2upper_passive_external_recv(blk,sub)
    type(blockdef), pointer :: blk
    procedure(extractarray) :: sub
end subroutine coor2upper_passive_external_recv

subroutine applyboundconds_passive_2d_corner_external(sub)
    procedure(extractarray) :: sub
    type(blockdef), pointer :: blk
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        select case (blk%loc_type)
        case (1)
            call xlyu_corner_passive_external(blk,sub)
        case (2)
        case (3)
            call xuyu_corner_passive_external(blk,sub)
        case (4)
        case (5)
        case (6)
        case (7)
            call xlyl_corner_passive_external(blk,sub)
        case (8)
        case (9)
            call xuyl_corner_passive_external(blk,sub)
        end select
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end subroutine applyboundconds_passive_2d_corner_external

subroutine xlyl_corner_passive_external(blk,sub)
    procedure(extractarray) :: sub
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
end subroutine xlyl_corner_passive_external

subroutine xuyl_corner_passive_external(blk,sub)
    procedure(extractarray) :: sub
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
end subroutine xuyl_corner_passive_external

subroutine xlyu_corner_passive_external(blk,sub)
    procedure(extractarray) :: sub
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
end subroutine xlyu_corner_passive_external

subroutine xuyu_corner_passive_external(blk,sub)
    procedure(extractarray) :: sub
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), pointer :: ap
end subroutine xuyu_corner_passive_external

end module boundary
