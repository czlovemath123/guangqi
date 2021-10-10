module communication
#include <petsc/finclude/petscsys.h>
use petscsys
use datastructure
use phylib
use eos
implicit none

contains

!the hash functions answer two questions, where from and where to

function hash_tag_send(blk_send,blk_recv,direction)
    !send tag has a simple hash function
    type(blockdef), pointer :: blk_send,blk_recv
    integer :: hash_tag_send,level_send,level_recv
    character(len=16) :: direction
    if (direction==east) then
        hash_tag_send=0
    else if (direction==south) then
        hash_tag_send=2
    else if (direction==west) then
        hash_tag_send=4
    else if (direction==north) then
        hash_tag_send=6
    else if (direction==se) then
        hash_tag_send=8
    else if (direction==sw) then
        hash_tag_send=9
    else if (direction==nw) then
        hash_tag_send=10
    else if (direction==ne) then
        hash_tag_send=11
    end if
    hash_tag_send=blk_send%blk_id*20+hash_tag_send
end function hash_tag_send

function hash_tag_recv(blk_send,blk_recv,direction)
    !need to consider different boundary conditions
    !all edges obey the same communication logic, corners are different
    !direction is with respect to the receiving block
    type(blockdef), pointer :: blk_send,blk_recv
    integer :: hash_tag_recv,level_send,level_recv,loc_type
    character(len=16) :: direction
    loc_type=blk_recv%loc_type
    if (direction==east) then
        hash_tag_recv=4
    else if (direction==south) then
        hash_tag_recv=6
    else if (direction==west) then
        hash_tag_recv=0
    else if (direction==north) then
        hash_tag_recv=2
    else if (direction==se) then
        if (loc_type==3.or.loc_type==6.or.loc_type==7.or.loc_type==8.or.loc_type==9) then
            if (loc_type==3.or.loc_type==6) then
                if (hydro_bound_type(2)==3) then
                    !periodic, not common
                    hash_tag_recv=10
                else
                    !transmissive or reflective
                    hash_tag_recv=11
                end if
            else if (loc_type==7.or.loc_type==8) then
                if (hydro_bound_type(3)==3) then
                    !periodic, cylindrical uses this boundary condition
                    hash_tag_recv=10
                else
                    !transmissive or reflective
                    hash_tag_recv=9
                end if
            else
                if (hydro_bound_type(2)==3.and.hydro_bound_type(3)==3) then
                    !both periodic, not common
                    hash_tag_recv=10
                else if (hydro_bound_type(2)==3) then
                    !not common
                    hash_tag_recv=9
                else if (hydro_bound_type(3)==3) then
                    !cylindrical
                    hash_tag_recv=11
                else
                    !no communication is needed
                end if
            end if
        else
            hash_tag_recv=10
        end if
    else if (direction==sw) then
        if (loc_type==1.or.loc_type==4.or.loc_type==7.or.loc_type==8.or.loc_type==9) then
            if (loc_type==1.or.loc_type==4) then
                if (hydro_bound_type(1)==3) then
                    !not common
                    hash_tag_recv=11
                else
                    hash_tag_recv=10
                end if
            else if (loc_type==8.or.loc_type==9) then
                if (hydro_bound_type(3)==3) then
                    !cylindrical
                    hash_tag_recv=11
                else
                    hash_tag_recv=8
                end if
            else
                if (hydro_bound_type(1)==3.and.hydro_bound_type(3)==3) then
                    !both periodic, not common
                    hash_tag_recv=11
                else if (hydro_bound_type(1)==3) then
                    !not common
                    hash_tag_recv=8
                else if (hydro_bound_type(3)==3) then
                    !cylindrical
                    hash_tag_recv=10
                else
                    !no communication is needed
                end if
            end if
        else
            hash_tag_recv=11
        end if
    else if (direction==nw) then
        if (loc_type==1.or.loc_type==2.or.loc_type==3.or.loc_type==4.or.loc_type==7) then
            if (loc_type==4.or.loc_type==7) then
                if (hydro_bound_type(1)==3) then
                    !not common
                    hash_tag_recv=8
                else
                    hash_tag_recv=9
                end if
            else if (loc_type==2.or.loc_type==3) then
                if (hydro_bound_type(4)==3) then
                    !cylindrical
                    hash_tag_recv=8
                else
                    hash_tag_recv=11
                end if
            else
                if (hydro_bound_type(1)==3.and.hydro_bound_type(4)==3) then
                    !both periodic, not common
                    hash_tag_recv=8
                else if (hydro_bound_type(1)==3) then
                    !not common
                    hash_tag_recv=11
                else if (hydro_bound_type(4)==3) then
                    !cylindrical
                    hash_tag_recv=9
                else
                    !no communication is needed
                end if
            end if
        else
            hash_tag_recv=8
        end if
    else if (direction==ne) then
        if (loc_type==1.or.loc_type==2.or.loc_type==3.or.loc_type==6.or.loc_type==9) then
            if (loc_type==1.or.loc_type==2) then
                if (hydro_bound_type(4)==3) then
                    !cylindrical
                    hash_tag_recv=9
                else
                    hash_tag_recv=10
                end if
            else if (loc_type==6.or.loc_type==9) then
                if (hydro_bound_type(2)==3) then
                    !not common
                    hash_tag_recv=9
                else
                    hash_tag_recv=8
                end if
            else
                if (hydro_bound_type(2)==3.and.hydro_bound_type(4)==3) then
                    !both periodic, not common
                    hash_tag_recv=9
                else if (hydro_bound_type(2)==3) then
                    !not common
                    hash_tag_recv=10
                else if (hydro_bound_type(4)==3) then
                    !cylindrical
                    hash_tag_recv=8
                else
                    !no communication is needed
                end if
            end if
        else
            hash_tag_recv=9
        end if
    end if
    hash_tag_recv=blk_send%blk_id*20+hash_tag_recv
end function hash_tag_recv

subroutine unlink_neighbour_blocks(blk)
    !remove the links between blocks
    type(blockdef), pointer :: blk
    type(blockneighbour) ,pointer :: blk_nb
    integer :: i
    if (nd==1) then
        if (associated(blk%nb_l)) then
            deallocate(blk%nb_l)
        end if
        if (associated(blk%nb_r)) then
            deallocate(blk%nb_r)
        end if
        if (associated(blk%blk_xl)) then
            call unlink_neighbour_blocks(blk%blk_xl)
        end if
        if (associated(blk%blk_xu)) then
            call unlink_neighbour_blocks(blk%blk_xu)
        end if
    else if (nd==2) then
    end if
end subroutine unlink_neighbour_blocks

!neighbour relation is critical in 2d and 3d problems, it lays the ground for communication

subroutine link_neighbour_blocks()
    !get the neighbour block cpu_rank and the global block key
    type(blockdef), pointer :: blk,blk_search1,blk_search2
    type(blockneighbour), pointer :: blk_nb
    integer :: i,key(3),key_nb1(3),key_nb2(3)
    if (nd==1) then
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            if (blk%loc_type==1) then
                call find_right_neighbour(blk)
            else if (blk%loc_type==2) then
                call find_left_neighbour(blk)
                call find_right_neighbour(blk)
            else if (blk%loc_type==3) then
                call find_left_neighbour(blk)
            end if
            if (i/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
        nullify(blk)
    else if (nd==2) then
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            call find_east_neighbour(blk)
            call find_south_neighbour(blk)
            call find_west_neighbour(blk)
            call find_north_neighbour(blk)
            call find_southeast_neighbour(blk)
            call find_southwest_neighbour(blk)
            call find_northwest_neighbour(blk)
            call find_northeast_neighbour(blk)
            if (i/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
        nullify(blk)
    else
    end if
end subroutine link_neighbour_blocks

subroutine examine_neighbour_relation()
end subroutine examine_neighbour_relation

subroutine establish_a_neighbourlink_1d(blk,direction)
    type(blockdef), pointer :: blk,blk_temp
    character(len=1) :: direction
    if (direction=='l') then
        blk_temp=>blk%blk_pre
        call new_neighbour_block(blk%nb_l)
        blk%nb_l%blk=>blk_temp
        call block_id_to_processor_rank(np_nblk,blk_temp%blk_id,blk%nb_l%cpu_rank)
    else if (direction=='r') then
        blk_temp=>blk%blk_next
        call new_neighbour_block(blk%nb_r)
        blk%nb_r%blk=>blk_temp
        call block_id_to_processor_rank(np_nblk,blk_temp%blk_id,blk%nb_r%cpu_rank)
    end if
end subroutine establish_a_neighbourlink_1d

subroutine find_left_neighbour(blk)
    type(blockdef), pointer :: blk
    character(len=1) :: direction
    direction='l'
    call establish_a_neighbourlink_1d(blk,direction)
end subroutine find_left_neighbour

subroutine find_right_neighbour(blk)
    type(blockdef), pointer :: blk
    character(len=1) :: direction
    direction='r'
    call establish_a_neighbourlink_1d(blk,direction)
end subroutine find_right_neighbour

subroutine link_a_neighbour_2d(blk,direction,blk_nb)
    !given a blk and a direction, return its neighbour on the direction
    !if no neighbouring block, return null
    type(blockdef), pointer :: blk
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    if (direction==east) then
        if (associated(blk%nb_e)) then
            blk_nb=>blk%nb_e
        else
            nullify(blk_nb)
        end if
    else if (direction==south) then
        if (associated(blk%nb_s)) then
            blk_nb=>blk%nb_s
        else
            nullify(blk_nb)
        end if
    else if (direction==west) then
        if (associated(blk%nb_w)) then
            blk_nb=>blk%nb_w
        else
            nullify(blk_nb)
        end if
    else if (direction==north) then
        if (associated(blk%nb_n)) then
            blk_nb=>blk%nb_n
        else
            nullify(blk_nb)
        end if
    else if (direction==se) then
        if (associated(blk%nb_se)) then
            blk_nb=>blk%nb_se
        else
            nullify(blk_nb)
        end if
    else if (direction==sw) then
        if (associated(blk%nb_sw)) then
            blk_nb=>blk%nb_sw
        else
            nullify(blk_nb)
        end if
    else if (direction==nw) then
        if (associated(blk%nb_nw)) then
            blk_nb=>blk%nb_nw
        else
            nullify(blk_nb)
        end if
    else if (direction==ne) then
        if (associated(blk%nb_ne)) then
            blk_nb=>blk%nb_ne
        else
            nullify(blk_nb)
        end if
    end if
end subroutine link_a_neighbour_2d

subroutine check_a_neighbour_2d(blk,direction,blk_nb,same)
    type(blockdef), pointer :: blk
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    logical :: same
    if (direction==east) then
        blk_nb=>blk%nb_e
    else if (direction==south) then
        blk_nb=>blk%nb_s
    else if (direction==west) then
        blk_nb=>blk%nb_w
    else if (direction==north) then
        blk_nb=>blk%nb_n
    else if (direction==se) then
        blk_nb=>blk%nb_se
    else if (direction==sw) then
        blk_nb=>blk%nb_sw
    else if (direction==nw) then
        blk_nb=>blk%nb_nw
    else if (direction==ne) then
        blk_nb=>blk%nb_ne
    end if
    if (blk_nb%cpu_rank==rank) then
        same=.true.
    else
        same=.false.
    end if
end subroutine check_a_neighbour_2d

subroutine establish_a_neighbourlink_2d(blk,key_nb,direction)
    !given blk, the key of its neighbour, and the relative direction of the neighbour
    !the inclusive subroutine that links a block to its neighbours
    type(blockdef), pointer :: blk,blk_search
    integer :: key(3),key_nb(3)
    character(len=16) :: direction
    key=blk%key
    call block_key_to_pointer(key_nb,blk_search)
    if (blk_search%key(3)-key(3)==-1) then
    else
        if (associated(blk_search%blk_xlyl)) then
        else
            if (direction==east) then
                call new_neighbour_block(blk%nb_e)
                blk%nb_e%blk=>blk_search
                call block_id_to_processor_rank(np_nblk,blk_search%blk_id,blk%nb_e%cpu_rank)
            else if (direction==south) then
                call new_neighbour_block(blk%nb_s)
                blk%nb_s%blk=>blk_search
                call block_id_to_processor_rank(np_nblk,blk_search%blk_id,blk%nb_s%cpu_rank)
            else if (direction==west) then
                call new_neighbour_block(blk%nb_w)
                blk%nb_w%blk=>blk_search
                call block_id_to_processor_rank(np_nblk,blk_search%blk_id,blk%nb_w%cpu_rank)
            else if (direction==north) then
                call new_neighbour_block(blk%nb_n)
                blk%nb_n%blk=>blk_search
                call block_id_to_processor_rank(np_nblk,blk_search%blk_id,blk%nb_n%cpu_rank)
            else if (direction==se) then
                call new_neighbour_block(blk%nb_se)
                blk%nb_se%blk=>blk_search
                call block_id_to_processor_rank(np_nblk,blk_search%blk_id,blk%nb_se%cpu_rank)
            else if (direction==sw) then
                call new_neighbour_block(blk%nb_sw)
                blk%nb_sw%blk=>blk_search
                call block_id_to_processor_rank(np_nblk,blk_search%blk_id,blk%nb_sw%cpu_rank)
            else if (direction==nw) then
                call new_neighbour_block(blk%nb_nw)
                blk%nb_nw%blk=>blk_search
                call block_id_to_processor_rank(np_nblk,blk_search%blk_id,blk%nb_nw%cpu_rank)
            else if (direction==ne) then
                call new_neighbour_block(blk%nb_ne)
                blk%nb_ne%blk=>blk_search
                call block_id_to_processor_rank(np_nblk,blk_search%blk_id,blk%nb_ne%cpu_rank)
            end if
        end if
    end if
end subroutine establish_a_neighbourlink_2d

subroutine find_east_neighbour(blk)
    !if there is an east neighbour, link blk%nb_e to the neighbour
    !else, nullify blk%nb_e
    type(blockdef), pointer :: blk,blk_search
    integer :: key(3),key_nb(3)
    character(len=16) :: direction
    key=blk%key
    direction=east
    if (blk%loc_type==3.or.blk%loc_type==6.or.blk%loc_type==9) then
        !the east neighbour is a part of the boundary condition
        if (hydro_bound_type(2)==3) then
            key_nb=(/1,key(2),key(3)/)
            call establish_a_neighbourlink_2d(blk,key_nb,direction)
        else
            nullify(blk%nb_e)
        end if
    else
        key_nb=(/key(1)+1,key(2),key(3)/)
        call establish_a_neighbourlink_2d(blk,key_nb,direction)
    end if
end subroutine find_east_neighbour

subroutine find_south_neighbour(blk)
    !if there is a south neighbour, link blk%nb_s to the neighbour
    !else, nullify blk%nb_s
    type(blockdef), pointer :: blk,blk_search
    integer :: key(3),key_nb(3)
    character(len=16) :: direction
    key=blk%key
    direction=south
    if (blk%loc_type==7.or.blk%loc_type==8.or.blk%loc_type==9) then
        !the south neighbour is a part of the boundary condition
        if (hydro_bound_type(3)==3) then
            key_nb=(/key(1),ny_blks*2**(key(3)-llevel_max),key(3)/)
            call establish_a_neighbourlink_2d(blk,key_nb,direction)
        else
            nullify(blk%nb_s)
        end if
    else
        key_nb=(/key(1),key(2)-1,key(3)/)
        call establish_a_neighbourlink_2d(blk,key_nb,direction)
    end if
end subroutine find_south_neighbour

subroutine find_west_neighbour(blk)
    !if there is a west neighbour, link blk%nb_w to the neighbour
    !else, nullify blk%nb_w
    type(blockdef), pointer :: blk,blk_search
    integer :: key(3),key_nb(3)
    character(len=16) :: direction
    key=blk%key
    direction=west
    if (blk%loc_type==1.or.blk%loc_type==4.or.blk%loc_type==7) then
        !the west neighbour is a part of the boundary condition
        if (hydro_bound_type(1)==3) then
            key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2),key(3)/)
            call establish_a_neighbourlink_2d(blk,key_nb,direction)
        else
            nullify(blk%nb_w)
        end if
    else
        key_nb=(/key(1)-1,key(2),key(3)/)
        call establish_a_neighbourlink_2d(blk,key_nb,direction)
    end if
end subroutine find_west_neighbour

subroutine find_north_neighbour(blk)
    !if there is a north neighbour, link blk%nb_n to the neighbour
    !else, nullify blk%nb_n
    type(blockdef), pointer :: blk,blk_search
    integer :: key(3),key_nb(3)
    character(len=16) :: direction
    key=blk%key
    direction=north
    if (blk%loc_type==1.or.blk%loc_type==2.or.blk%loc_type==3) then
        !the north neighbour is a part of the boundary condition
        if (hydro_bound_type(4)==3) then
            key_nb=(/key(1),1,key(3)/)
            call establish_a_neighbourlink_2d(blk,key_nb,direction)
        else
            nullify(blk%nb_n)
        end if
    else
        key_nb=(/key(1),key(2)+1,key(3)/)
        call establish_a_neighbourlink_2d(blk,key_nb,direction)
    end if
end subroutine find_north_neighbour

subroutine find_southeast_neighbour(blk)
    type(blockdef), pointer :: blk,blk_search
    integer :: key(3),key_nb(3),loc_type
    character(len=16) :: direction
    key=blk%key
    direction=se
    loc_type=blk%loc_type
    if (loc_type==3.or.loc_type==6.or.loc_type==7.or.loc_type==8.or.loc_type==9) then
        if (loc_type==3.or.loc_type==6) then
            if (hydro_bound_type(2)==3) then
                key_nb=(/1,key(2)-1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(2)==9) then
                nullify(blk%nb_se)
            else
                key_nb=(/key(1),key(2)-1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            end if
        end if
        if (loc_type==7.or.loc_type==8) then
            if (hydro_bound_type(3)==3) then
                key_nb=(/key(1)+1,ny_blks*2**(key(3)-llevel_max),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(3)==9) then
                nullify(blk%nb_se)
            else
                key_nb=(/key(1)+1,key(2),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            end if
        end if
        if (loc_type==9) then
            if (hydro_bound_type(2)==9.or.hydro_bound_type(3)==9) then
                nullify(blk%nb_se)
            else if (hydro_bound_type(2)==3.and.hydro_bound_type(3)==3) then
            else if (hydro_bound_type(2)==3) then
                key_nb=(/1,key(2),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(3)==3) then
                key_nb=(/key(1),ny_blks*2**(key(3)-llevel_max),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else
                nullify(blk%nb_se)
            end if
        end if
    else
        key_nb=(/key(1)+1,key(2)-1,key(3)/)
        call establish_a_neighbourlink_2d(blk,key_nb,direction)
    end if
end subroutine find_southeast_neighbour

subroutine find_southwest_neighbour(blk)
    type(blockdef), pointer :: blk,blk_search
    integer :: key(3),key_nb(3),loc_type
    character(len=16) :: direction
    key=blk%key
    direction=sw
    loc_type=blk%loc_type
    if (loc_type==1.or.loc_type==4.or.loc_type==7.or.loc_type==8.or.loc_type==9) then
        if (loc_type==1.or.loc_type==4) then
            if (hydro_bound_type(1)==3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2)-1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(1)==9) then
                nullify(blk%nb_sw)
            else
                key_nb=(/key(1),key(2)-1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            end if
        end if
        if (loc_type==8.or.loc_type==9) then
            if (hydro_bound_type(3)==3) then
                key_nb=(/key(1)-1,ny_blks*2**(key(3)-llevel_max),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(3)==9) then
                nullify(blk%nb_sw)
            else
                key_nb=(/key(1)-1,key(2),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            end if
        end if
        if (loc_type==7) then
            if (hydro_bound_type(1)==9.or.hydro_bound_type(3)==9) then
                nullify(blk%nb_sw)
            else if (hydro_bound_type(1)==3.and.hydro_bound_type(3)==3) then
            else if (hydro_bound_type(1)==3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(3)==3) then
                key_nb=(/key(1),ny_blks*2**(key(3)-llevel_max),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else
                nullify(blk%nb_sw)
            end if
        end if
    else
        key_nb=(/key(1)-1,key(2)-1,key(3)/)
        call establish_a_neighbourlink_2d(blk,key_nb,direction)
    end if
end subroutine find_southwest_neighbour

subroutine find_northwest_neighbour(blk)
    type(blockdef), pointer :: blk,blk_search
    integer :: key(3),key_nb(3),loc_type
    character(len=16) :: direction
    key=blk%key
    direction=nw
    loc_type=blk%loc_type
    if (loc_type==1.or.loc_type==2.or.loc_type==3.or.loc_type==4.or.loc_type==7) then
        if (loc_type==4.or.loc_type==7) then
            if (hydro_bound_type(1)==3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2)+1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(1)==9) then
                nullify(blk%nb_nw)
            else
                key_nb=(/key(1),key(2)+1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            end if
        end if
        if (loc_type==2.or.loc_type==3) then
            if (hydro_bound_type(4)==3) then
                key_nb=(/key(1)-1,1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(4)==9) then
                nullify(blk%nb_nw)
            else
                key_nb=(/key(1)-1,key(2),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            end if
        end if
        if (loc_type==1) then
            if (hydro_bound_type(1)==9.and.hydro_bound_type(4)==9) then
            else if (hydro_bound_type(1)==3.and.hydro_bound_type(4)==3) then
            else if (hydro_bound_type(1)==3) then
                key_nb=(/nx_blks*2**(key(3)-llevel_max),key(2),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(4)==3) then
                key_nb=(/key(1),1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else
                nullify(blk%nb_nw)
            end if
        end if
    else
        key_nb=(/key(1)-1,key(2)+1,key(3)/)
        call establish_a_neighbourlink_2d(blk,key_nb,direction)
    end if
end subroutine find_northwest_neighbour

subroutine find_northeast_neighbour(blk)
    type(blockdef), pointer :: blk,blk_search
    integer :: key(3),key_nb(3),loc_type
    character(len=16) :: direction
    key=blk%key
    direction=ne
    loc_type=blk%loc_type
    if (loc_type==1.or.loc_type==2.or.loc_type==3.or.loc_type==6.or.loc_type==9) then
        if (loc_type==6.or.loc_type==9) then
            if (hydro_bound_type(2)==3) then
                key_nb=(/1,key(2)+1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(2)==9) then
                nullify(blk%nb_ne)
            else
                key_nb=(/key(1),key(2)+1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            end if
        end if
        if (loc_type==1.or.loc_type==2) then
            if (hydro_bound_type(4)==3) then
                key_nb=(/key(1)+1,1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(4)==9) then
                nullify(blk%nb_ne)
            else
                key_nb=(/key(1)+1,key(2),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            end if
        end if
        if (loc_type==3) then
            if (hydro_bound_type(2)==9.and.hydro_bound_type(4)==9) then
            else if (hydro_bound_type(2)==3.and.hydro_bound_type(4)==3) then
            else if (hydro_bound_type(2)==3) then
                key_nb=(/1,key(2),key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else if (hydro_bound_type(4)==3) then
                key_nb=(/key(1),1,key(3)/)
                call establish_a_neighbourlink_2d(blk,key_nb,direction)
            else
                nullify(blk%nb_ne)
            end if
        end if
    else
        key_nb=(/key(1)+1,key(2)+1,key(3)/)
        call establish_a_neighbourlink_2d(blk,key_nb,direction)
    end if
end subroutine find_northeast_neighbour

subroutine corner_mirror_reflect_hydro(a,mirror,reflect,b)
    !the corner cells may need mirror and reflect operation at the boundaries
    real(8), dimension(2,2,5) :: a,b
    character(len=1) :: mirror
    logical :: reflect
    integer :: i,j
    if (mirror=='x') then
        if (reflect) then
            do j=1,2
                do i=1,2
                    b(i,j,1:5)=a(3-i,j,1:5)
                    b(i,j,2)=-a(3-i,j,2)
                end do
            end do
        else
            do j=1,2
                do i=1,2
                    b(i,j,1:5)=a(3-i,j,1:5)
                end do
            end do
        end if
    else if (mirror=='y') then
        if (reflect) then
            do j=1,2
                do i=1,2
                    b(i,j,1:5)=a(i,3-j,1:5)
                    b(i,j,3)=-a(i,3-j,3)
                end do
            end do
        else
            do j=1,2
                do i=1,2
                    b(i,j,1:5)=a(i,3-j,1:5)
                end do
            end do
        end if
    end if
end subroutine corner_mirror_reflect_hydro

subroutine corner_mirror_scalar(a,mirror,b)
    !scalar does not have reflect operation
    real(8), dimension(2,2) :: a,b
    character(len=1) :: mirror
    integer :: i,j
    if (mirror=='x') then
        do j=1,2
            do i=1,2
                b(i,j)=a(3-i,j)
            end do
        end do
    else if (mirror=='y') then
        do j=1,2
            do i=1,2
                b(i,j)=a(i,3-j)
            end do
        end do
    end if
end subroutine corner_mirror_scalar

subroutine communicate_flux()
    !inter core communication
    type(blockdef), pointer :: blk_traversal,blk_xl,blk_xu
    real(8) :: msg_out(5),msg_in(5),xflux_xl(5),xflux_xu(5)
    integer :: level,level_pre,level_next,i,tag,ierr,stat(MPI_STATUS_SIZE),req
    if (nd==1) then
        xflux_xl=0d0;xflux_xu=0d0
        if (np>1) then
            if (rank==0) then
                msg_out(1:5)=blk_processor_tail%xflux(blk_size_nx,1,1,1:5)
                call mpi_isend(msg_out,5,MPI_REAL8,rank+1,1,MPI_COMM_WORLD,req,ierr)
            else if (rank==np-1) then
                msg_out(1:5)=blk_processor_head%xflux(0,1,1,1:5)
                call mpi_isend(msg_out,5,MPI_REAL8,rank-1,1,MPI_COMM_WORLD,req,ierr)
            else
                msg_out(1:5)=blk_processor_tail%xflux(blk_size_nx,1,1,1:5)
                call mpi_isend(msg_out,5,MPI_REAL8,rank+1,1,MPI_COMM_WORLD,req,ierr)
                msg_out(1:5)=blk_processor_head%xflux(0,1,1,1:5)
                call mpi_isend(msg_out,5,MPI_REAL8,rank-1,1,MPI_COMM_WORLD,req,ierr)
            end if
            if (rank==0) then
                call mpi_recv(xflux_xu,5,MPI_REAL8,rank+1,1,MPI_COMM_WORLD,stat,ierr)
            else if (rank==np-1) then
                call mpi_recv(xflux_xl,5,MPI_REAL8,rank-1,1,MPI_COMM_WORLD,stat,ierr)
            else
                call mpi_recv(xflux_xu,5,MPI_REAL8,rank+1,1,MPI_COMM_WORLD,stat,ierr)
                call mpi_recv(xflux_xl,5,MPI_REAL8,rank-1,1,MPI_COMM_WORLD,stat,ierr)
            end if
        end if
        blk_traversal=>blk_processor_head
        call communicate_xflux_1d_xl(blk_traversal,xflux_xl)
        blk_traversal=>blk_traversal%blk_next
        do i=2,np_nblk(rank+1)-1
            call communicate_xflux_1d_interior(blk_traversal)
            blk_traversal=>blk_traversal%blk_next
        end do
        call communicate_xflux_1d_xu(blk_traversal,xflux_xu)
        nullify(blk_traversal)
    else if (nd==2) then
    end if
end subroutine communicate_flux

subroutine communicate_xflux_1d_xl(blk,xflux_xl)
    type(blockdef), pointer :: blk,blk_xu
    real(8) :: xflux_xl(5)
    integer :: level,level_pre,level_next
    level=blk%level
    if (rank==0) then
        blk_xu=>blk%blk_next
        level_next=blk_xu%level
        if (level<level_next) then
            blk%xflux(blk_size_nx,1,1,1:5)=blk_xu%xflux(0,1,1,1:5)
        end if
    else
        level_pre=blk%blk_pre%level
        blk_xu=>blk%blk_next
        level_next=blk_xu%level
        if (level<level_pre) then
            blk%xflux(0,1,1,1:5)=xflux_xl
        end if
        if (level<level_next) then
            blk%xflux(blk_size_nx,1,1,1:5)=blk_xu%xflux(0,1,1,1:5)
        end if
    end if
end subroutine communicate_xflux_1d_xl

subroutine communicate_xflux_1d_interior(blk)
    type(blockdef), pointer :: blk,blk_xu
    integer :: level,level_next
    level=blk%level
    blk_xu=>blk%blk_next
    level_next=blk_xu%level
    !the interior blocks only need to check every interface to its right
    if (level>level_next) then
        blk_xu%xflux(0,1,1,1:5)=blk%xflux(blk_size_nx,1,1,1:5)
    else if (level<level_next) then
        blk%xflux(blk_size_nx,1,1,1:5)=blk_xu%xflux(0,1,1,1:5)
    end if
end subroutine communicate_xflux_1d_interior

subroutine communicate_xflux_1d_xu(blk,xflux_xu)
    type(blockdef), pointer :: blk,blk_xl
    real(8) :: xflux_xu(5)
    integer :: level,level_pre,level_next
    level=blk%level
    if (rank==np-1) then
        blk_xl=>blk%blk_pre
        level_pre=blk_xl%level
        if (level<level_pre) then
            blk%xflux(0,1,1,1:5)=blk_xl%xflux(blk_size_nx,1,1,1:5)
        end if
    else
        blk_xl=>blk%blk_pre
        level_pre=blk%blk_pre%level
        level_next=blk%blk_next%level
        if (level<level_pre) then
            blk%xflux(0,1,1,1:5)=blk_xl%xflux(blk_size_nx,1,1,1:5)
        end if
        if (level<level_next) then
            blk%xflux(blk_size_nx,1,1,1:5)=xflux_xu
        end if
    end if
end subroutine communicate_xflux_1d_xu

subroutine boundary_restriction_1d(blk1,ibound,w,u,egv,temp)
    !do restriction on one boundary size of blk1, the result is stored in w,u,egv,temp
    type(blockdef), pointer, intent(in) :: blk1
    real(8) :: w(5,2),u(5,2),temp(2),egv(2),vol1,vol2
    integer :: ibound,i,j
    u=0d0;w=0d0;temp=0d0;egv=0d0
    select case(ibound)
    case(1)         !xl boundary
        do i=1,n_hydro_guard
            if (igeometry==0) then
                u(1:5,i)=(blk1%u(2*i-1,1,1,1:5)+blk1%u(2*i,1,1,1:5))/2
            else if (igeometry==2) then
                u(1:5,i)=(blk1%u(2*i-1,1,1,1:5)*blk1%vol(2*i-1,1,1)+blk1%u(2*i,1,1,1:5)*blk1%vol(2*i,1,1))/(blk1%vol(2*i-1,1,1)+blk1%vol(2*i,1,1))
            end if
#if         ieos==1
            call utow(u(1:5,i),w(1:5,i))
            temp(i)=calculate_temp_from_p_rho(w(5,i),w(1,i))
            egv(i)=w(5,i)/(gamma_gas-one)
#elif       ieos==2
            call eos_hllc_analytic_utow(u(1:5,i),w(1:5,i),temp(i),egv(i))
#endif
        end do
    case(2)         !xu boundary
        do i=1,n_hydro_guard
            if (igeometry==0) then
                u(1:5,i)=(blk1%u(blk_size_nx-2*n_hydro_guard+2*i-1,1,1,1:5)   &
                    +blk1%u(blk_size_nx-2*n_hydro_guard+2*i,1,1,1:5))/2
            else if (igeometry==2) then
                vol1=blk1%vol(blk_size_nx-2*n_hydro_guard+2*i-1,1,1)
                vol2=blk1%vol(blk_size_nx-2*n_hydro_guard+2*i,1,1)
                u(1:5,i)=(blk1%u(blk_size_nx-2*n_hydro_guard+2*i-1,1,1,1:5)*vol1   &
                    +blk1%u(blk_size_nx-2*n_hydro_guard+2*i,1,1,1:5)*vol2)/(vol1+vol2)
            end if
#if         ieos==1
            call utow(u(1:5,i),w(1:5,i))
            temp(i)=calculate_temp_from_p_rho(w(5,i),w(1,i))
            egv(i)=w(5,i)/(gamma_gas-one)
#elif       ieos==2
            call eos_hllc_analytic_utow(u(1:5,i),w(1:5,i),temp(i),egv(i))
#endif
        end do
    case(3)         !yl boundary
    case(4)         !yu boundary
    end select
end subroutine boundary_restriction_1d

subroutine boundary_prolongation_1d(blk1,ibound,w,u,egv,temp)
    !do prolongation on one boundary size of blk1, the result is stored in w,u,egv,temp
    type(blockdef), pointer, intent(in) :: blk1
    real(8) :: w(5,2),u(5,2),temp(2),egv(2)
    real(8) :: diff_l,diff_r,diff_c,diff
    integer :: ibound,i,j
    u=0d0;w=0d0;temp=0d0;egv=0d0
    select case(ibound)
    case(1)         !xl boundary
        if (n_hydro_guard/=2) then
            print *,'not implemented'
            stop
        end if
        do j=1,5
            diff_l=blk1%w(1,1,1,j)-blk1%w(0,1,1,j)
            diff_r=blk1%w(2,1,1,j)-blk1%w(1,1,1,j)
            diff_c=(blk1%w(2,1,1,j)-blk1%w(0,1,1,j))/2
            if (diff_c/=zero) then
                diff=abs(diff_c)/diff_c*min(abs(diff_l),abs(diff_r),abs(diff_c))
            else
                diff=0d0
            end if
            diff=0d0
            w(j,1)=blk1%w(1,1,1,j)-half*diff
            w(j,2)=blk1%w(1,1,1,j)+half*diff
        end do
        do i=1,2
#if         ieos==1
            call wtou(w(1:5,i),u(1:5,i))
            temp(i)=calculate_temp_from_p_rho(w(5,i),w(1,i))
            egv(i)=w(5,i)/(gamma_gas-one)
#elif       ieos==2
            temp(i)=solvetp(w(5,i),w(1,i))
            egv(i)=egvrhot(w(1,i),temp(i))
            call eos_hllc_analytic_wtou(w(1:5,i),egv(i),u(1:5,i))
#endif
        end do
    case(2)         !xu boundary
        if (n_hydro_guard/=2) then
            print *,'not implemented'
            stop
        end if
        do j=1,5
            diff_l=blk1%w(blk_size_nx,1,1,j)-blk1%w(blk_size_nx-1,1,1,j)
            diff_r=blk1%w(blk_size_nx+1,1,1,j)-blk1%w(blk_size_nx,1,1,j)
            diff_c=(blk1%w(blk_size_nx+1,1,1,j)-blk1%w(blk_size_nx-1,1,1,j))/2
            if (diff_c/=zero) then
                diff=abs(diff_c)/diff_c*min(abs(diff_l),abs(diff_r),abs(diff_c))
            else
                diff=0d0
            end if
            diff=0d0
            w(j,1)=blk1%w(blk_size_nx,1,1,j)-half*diff
            w(j,2)=blk1%w(blk_size_nx,1,1,j)+half*diff
        end do
        do i=1,2
#if         ieos==1
            call wtou(w(1:5,i),u(1:5,i))
            temp(i)=calculate_temp_from_p_rho(w(5,i),w(1,i))
            egv(i)=w(5,i)/(gamma_gas-one)
#elif       ieos==2
            temp(i)=solvetp(w(5,i),w(1,i))
            egv(i)=egvrhot(w(1,i),temp(i))
            call eos_hllc_analytic_wtou(w(1:5,i),egv(i),u(1:5,i))
#endif
        end do
    case(3)         !yl boundary
    case(4)         !yu boundary
    end select
end subroutine boundary_prolongation_1d

subroutine communicate_hydro()
    !first restriction, then prolongation
    !communicate non-boundary edges and vertexes
    type(blockdef), pointer :: blk_traversal,blk_xl,blk_xu
    integer :: i,ierr
    if (nd==1) then
        call communicate_blocks_hydro_internal_1d()
        call communicate_blocks_hydro_intercores_1d()
    else if (nd==2) then
        call communicate_hydro_internal_2d()
        call communicate_hydro_external_2d()
    end if
end subroutine communicate_hydro

subroutine communicate_blocks_hydro_internal_1d()
    type(blockdef), pointer :: blk
    integer :: i,ierr
    character(len=128) :: str
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        select case (blk%loc_type)
        case(1)
            call communicate_hydro_right_internal(blk)
        case(2)
            call communicate_hydro_left_internal(blk)
            call communicate_hydro_right_internal(blk)
        case(3)
            call communicate_hydro_left_internal(blk)
        case default
            str='no such loc_type in 1D'
            print *,blk%blk_id,blk%loc_type
            call abort_achilles(str)
        end select
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end subroutine communicate_blocks_hydro_internal_1d

subroutine communicate_hydro_left_internal(blk)
    !restrict or prolongate or copy the left neighbour's boundary cells and fill blk's guard cells
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    real(8) :: w(5,2),u(5,2),temp(2),egv(2)
    integer :: i,nb_level,level,neighbour_rank
    if (associated(blk%nb_l)) then
        blk_nb=>blk%nb_l
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank==rank) then
            blk_temp=>blk_nb%blk
            level=blk%level
            nb_level=blk_temp%level
            if (nb_level>level) then
                call boundary_restriction_1d(blk_temp,2,w,u,egv,temp)
            else if (nb_level<level) then
                call boundary_prolongation_1d(blk_temp,2,w,u,egv,temp)
            else
                !same level
                w=transpose(blk_temp%w(blk_size_nx-1:blk_size_nx,1,1,1:5))
                u=transpose(blk_temp%u(blk_size_nx-1:blk_size_nx,1,1,1:5))
                egv=blk_temp%egv(blk_size_nx-1:blk_size_nx,1,1)
                temp=blk_temp%temp(blk_size_nx-1:blk_size_nx,1,1)
            end if
            blk%w(-1:0,1,1,1:5)=transpose(w)
            blk%u(-1:0,1,1,1:5)=transpose(u)
            blk%temp(-1:0,1,1)=temp
            blk%egv(-1:0,1,1)=egv
        end if
        nullify(blk_nb,blk_temp)
    end if
end subroutine communicate_hydro_left_internal

subroutine communicate_hydro_right_internal(blk)
    !restrict or prolongate or copy the right neighbour's boundary cells and fill blk's guard cells
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    real(8) :: w(5,2),u(5,2),temp(2),egv(2)
    integer :: i,nb_level,level,neighbour_rank
    if (associated(blk%nb_r)) then
        blk_nb=>blk%nb_r
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank==rank) then
            blk_temp=>blk_nb%blk
            level=blk%level
            nb_level=blk_temp%level
            if (nb_level>level) then
                call boundary_restriction_1d(blk_temp,1,w,u,egv,temp)
            else if (nb_level<level) then
                call boundary_prolongation_1d(blk_temp,1,w,u,egv,temp)
            else
                !same level
                w=transpose(blk_temp%w(1:2,1,1,1:5))
                u=transpose(blk_temp%u(1:2,1,1,1:5))
                egv=blk_temp%egv(1:2,1,1)
                temp=blk_temp%temp(1:2,1,1)
            end if
            blk%w(blk_size_nx+1:blk_size_nx+2,1,1,1:5)=transpose(w)
            blk%u(blk_size_nx+1:blk_size_nx+2,1,1,1:5)=transpose(u)
            blk%temp(blk_size_nx+1:blk_size_nx+2,1,1)=temp
            blk%egv(blk_size_nx+1:blk_size_nx+2,1,1)=egv
        end if
        nullify(blk_nb,blk_temp)
    end if
end subroutine communicate_hydro_right_internal

subroutine communicate_blocks_hydro_intercores_1d()
    type(blockdef), pointer :: blk
    integer :: i,ierr
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        select case (blk%loc_type)
        case(1)
            call communicate_right_neighbour_intercores(blk)
        case(2)
            call communicate_left_neighbour_intercores(blk)
            call communicate_right_neighbour_intercores(blk)
        case(3)
            call communicate_left_neighbour_intercores(blk)
        end select
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
end subroutine communicate_blocks_hydro_intercores_1d

subroutine communicate_left_neighbour_intercores(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    real(8) :: w(5,2),u(5,2),temp(2),egv(2),msg_out(12,2),msg_in(12,2),msg_tp(2,12)
    integer :: i,ierr,nb_level,level,neighbour_rank,send_count,recv_count,tag,req,stat(MPI_STATUS_SIZE)
    if (associated(blk%nb_l)) then
        blk_nb=>blk%nb_l
        blk_temp=>blk_nb%blk
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            level=blk%level
            nb_level=blk_temp%level
            if (nb_level>level) then
                call boundary_prolongation_1d(blk,1,w,u,egv,temp)
                msg_out(1:5,1:2)=w
                msg_out(6:10,1:2)=u
                msg_out(11,1:2)=egv
                msg_out(12,1:2)=temp
            else if (nb_level<level) then
                call boundary_restriction_1d(blk,1,w,u,egv,temp)
                msg_out(1:5,1:2)=w
                msg_out(6:10,1:2)=u
                msg_out(11,1:2)=egv
                msg_out(12,1:2)=temp
            else
                !same level
                msg_out(1:5,1:2)=transpose(blk%w(1:2,1,1,1:5))
                msg_out(6:10,1:2)=transpose(blk%u(1:2,1,1,1:5))
                msg_out(11,1:2)=blk%egv(1:2,1,1)
                msg_out(12,1:2)=blk%temp(1:2,1,1)
            end if
            send_count=size(msg_out)
            tag=1;call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
            call mpi_recv(msg_in,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
            msg_tp=transpose(msg_in)
            blk%w(-1:0,1,1,1:5)=msg_tp(1:2,1:5)
            blk%u(-1:0,1,1,1:5)=msg_tp(1:2,6:10)
            blk%egv(-1:0,1,1)=msg_tp(1:2,11)
            blk%temp(-1:0,1,1)=msg_tp(1:2,12)
        end if
    end if
end subroutine communicate_left_neighbour_intercores

subroutine communicate_right_neighbour_intercores(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    real(8) :: w(5,2),u(5,2),temp(2),egv(2),msg_out(12,2),msg_in(12,2),msg_tp(2,12)
    integer :: i,ierr,nb_level,level,neighbour_rank,send_count,recv_count,tag,req,stat(MPI_STATUS_SIZE)
    if (associated(blk%nb_r)) then
        blk_nb=>blk%nb_r
        blk_temp=>blk_nb%blk
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            level=blk%level
            nb_level=blk_temp%level
            if (nb_level>level) then
                call boundary_prolongation_1d(blk,2,w,u,egv,temp)
                msg_out(1:5,1:2)=w
                msg_out(6:10,1:2)=u
                msg_out(11,1:2)=egv
                msg_out(12,1:2)=temp
            else if (nb_level<level) then
                call boundary_restriction_1d(blk,2,w,u,egv,temp)
                msg_out(1:5,1:2)=w
                msg_out(6:10,1:2)=u
                msg_out(11,1:2)=egv
                msg_out(12,1:2)=temp
            else
                !same level
                msg_out(1:5,1:2)=transpose(blk%w(blk_size_nx-1:blk_size_nx,1,1,1:5))
                msg_out(6:10,1:2)=transpose(blk%u(blk_size_nx-1:blk_size_nx,1,1,1:5))
                msg_out(11,1:2)=blk%egv(blk_size_nx-1:blk_size_nx,1,1)
                msg_out(12,1:2)=blk%temp(blk_size_nx-1:blk_size_nx,1,1)
            end if
            send_count=size(msg_out)
            tag=1;call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
            call mpi_recv(msg_in,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
            msg_tp=transpose(msg_in)
            blk%w(blk_size_nx+1:blk_size_nx+2,1,1,1:5)=msg_tp(1:2,1:5)
            blk%u(blk_size_nx+1:blk_size_nx+2,1,1,1:5)=msg_tp(1:2,6:10)
            blk%egv(blk_size_nx+1:blk_size_nx+2,1,1)=msg_tp(1:2,11)
            blk%temp(blk_size_nx+1:blk_size_nx+2,1,1)=msg_tp(1:2,12)
        end if
    end if
end subroutine communicate_right_neighbour_intercores

subroutine communicate_hydro_internal_2d()
    type(blockdef), pointer :: blk
    integer :: i,ierr
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        call communicate_hydro_edge_internal(blk,east)
        call communicate_hydro_edge_internal(blk,south)
        call communicate_hydro_edge_internal(blk,west)
        call communicate_hydro_edge_internal(blk,north)
        call communicate_hydro_corner_internal(blk,se)
        call communicate_hydro_corner_internal(blk,sw)
        call communicate_hydro_corner_internal(blk,nw)
        call communicate_hydro_corner_internal(blk,ne)
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    nullify(blk)
end subroutine communicate_hydro_internal_2d

subroutine communicate_hydro_edge_internal(blk,direction)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    integer :: i,j,neighbour_rank
    call link_a_neighbour_2d(blk,direction,blk_nb)
    if (associated(blk_nb)) then
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank==rank) then
            blk_temp=>blk_nb%blk
            if (direction==east) then
                blk%w(blk_size_nx+1:blk_size_nx+2,1:blk_size_ny,1,1:5)=blk_temp%w(1:2,1:blk_size_ny,1,1:5)
            else if (direction==south) then
                blk%w(1:blk_size_nx,-1:0,1,1:5)=blk_temp%w(1:blk_size_nx,blk_size_ny-1:blk_size_ny,1,1:5)
            else if (direction==west) then
                blk%w(-1:0,1:blk_size_ny,1,1:5)=blk_temp%w(blk_size_nx-1:blk_size_nx,1:blk_size_ny,1,1:5)
            else if (direction==north) then
                blk%w(1:blk_size_nx,blk_size_ny+1:blk_size_ny+2,1,1:5)=blk_temp%w(1:blk_size_nx,1:2,1,1:5)
            end if
            nullify(blk_temp)
        end if
    end if
end subroutine communicate_hydro_edge_internal

subroutine extract_corner_hydro(blk,direction,w)
    type(blockdef), pointer :: blk
    character(len=16) :: direction
    real(8), dimension(2,2,5) :: w
    if (direction==se) then
        w=blk%w(blk_size_nx-1:blk_size_nx,1:2,1,1:5)
    else if (direction==sw) then
        w=blk%w(1:2,1:2,1,1:5)
    else if (direction==nw) then
        w=blk%w(1:2,blk_size_ny-1:blk_size_ny,1,1:5)
    else if (direction==ne) then
        w=blk%w(blk_size_nx-1:blk_size_nx,blk_size_ny-1:blk_size_ny,1,1:5)
    end if
end subroutine extract_corner_hydro

subroutine communicate_hydro_corner_internal(blk,direction)
    !only operate on corners that need communication
    !for user defined boundary condition and domain corners, see corner boundary conditions
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    integer :: i,j,neighbour_rank
    real(8), dimension(2,2,5) :: w_in,u_in,w_out,u_out
    real(8), dimension(2,2) :: temp_in,egv_in,temp_out,egv_out
    character(len=16) :: direction
    character(len=1) :: mirror
    logical :: reflect
    call link_a_neighbour_2d(blk,direction,blk_nb)
    if (associated(blk_nb)) then
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank==rank) then
            blk_temp=>blk_nb%blk
            select case (blk%loc_type)
            case (1)
                !corner block, transmissive, reflective, and periodic
                if (direction==ne) then
                    if (hydro_bound_type(4)==3) then
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==se) then
                    call extract_corner_hydro(blk_temp,nw,w_in)
                    w_out=w_in
                else if (direction==sw) then
                    if (hydro_bound_type(1)==3) then
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==nw) then
                    if (hydro_bound_type(1)==3.and.hydro_bound_type(4)/=3) then
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(1)/=3.and.hydro_bound_type(4)==3) then
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(1)==3.and.hydro_bound_type(4)==3) then
                        call extract_corner_hydro(blk_temp,se,w_in)
                        w_out=w_in
                    end if
                end if
            case (2)
                !corner on edge, transmissive, reflective, periodic or user specified
                if (direction==ne) then
                    if (hydro_bound_type(4)==3) then
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==se) then
                    call extract_corner_hydro(blk_temp,nw,w_in)
                    w_out=w_in
                else if (direction==sw) then
                    call extract_corner_hydro(blk_temp,ne,w_in)
                    w_out=w_in
                else if (direction==nw) then
                    if (hydro_bound_type(4)==3) then
                        call extract_corner_hydro(blk_temp,se,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                end if
            case (3)
                !corner block, transmissive, reflective, periodic, and user specified
                if (direction==ne) then
                    if (hydro_bound_type(2)==3.and.hydro_bound_type(4)/=3) then
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(2)/=3.and.hydro_bound_type(4)==3) then
                        call extract_corner_hydro(blk_temp,se,w_in)
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(2)==3.and.hydro_bound_type(4)==3) then
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        w_out=w_in
                    end if
                else if (direction==se) then
                    if (hydro_bound_type(2)==3) then
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==sw) then
                    call extract_corner_hydro(blk_temp,ne,w_in)
                    w_out=w_in
                else if (direction==nw) then
                    if (hydro_bound_type(4)==3) then
                        call extract_corner_hydro(blk_temp,se,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                end if
            case (4)
                !corner on edge, transmissive, reflective, periodic or user specified
                if (direction==ne) then
                    call extract_corner_hydro(blk_temp,sw,w_in)
                    w_out=w_in
                else if (direction==se) then
                    call extract_corner_hydro(blk_temp,nw,w_in)
                    w_out=w_in
                else if (direction==sw) then
                    if (hydro_bound_type(1)==3) then
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==nw) then
                    if (hydro_bound_type(1)==3) then
                        call extract_corner_hydro(blk_temp,se,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                end if
            case (5)
                !inside
                if (direction==ne) then
                    call extract_corner_hydro(blk_temp,sw,w_in)
                else if (direction==se) then
                    call extract_corner_hydro(blk_temp,nw,w_in)
                else if (direction==sw) then
                    call extract_corner_hydro(blk_temp,ne,w_in)
                else if (direction==nw) then
                    call extract_corner_hydro(blk_temp,se,w_in)
                end if
                w_out=w_in
            case (6)
                !corner on edge, transmissive, reflective, periodic or user specified
                if (direction==ne) then
                    if (hydro_bound_type(2)==3) then
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,se,w_in)
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==se) then
                    if (hydro_bound_type(2)==3) then
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==sw) then
                    call extract_corner_hydro(blk_temp,ne,w_in)
                    w_out=w_in
                else if (direction==nw) then
                    call extract_corner_hydro(blk_temp,se,w_in)
                    w_out=w_in
                end if
            case (7)
                !corner block, transmissive, reflective, periodic, and user specified
                if (direction==ne) then
                    call extract_corner_hydro(blk_temp,sw,w_in)
                    w_out=w_in
                else if (direction==se) then
                    if (hydro_bound_type(3)==3) then
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==sw) then
                    if (hydro_bound_type(1)==3.and.hydro_bound_type(3)/=3) then
                        call extract_corner_hydro(blk_temp,se,w_in)
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(1)/=3.and.hydro_bound_type(3)==3) then
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(1)==3.and.hydro_bound_type(3)==3) then
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        w_out=w_in
                    end if
                else if (direction==nw) then
                    if (hydro_bound_type(1)==3) then
                        call extract_corner_hydro(blk_temp,se,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                end if
            case (8)
                !corner on edge, transmissive, reflective, periodic or user specified
                if (direction==ne) then
                    call extract_corner_hydro(blk_temp,sw,w_in)
                    w_out=w_in
                else if (direction==se) then
                    if (hydro_bound_type(3)==3) then
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==sw) then
                    if (hydro_bound_type(3)==3) then
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,se,w_in)
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==nw) then
                    call extract_corner_hydro(blk_temp,se,w_in)
                    w_out=w_in
                end if
            case (9)
                !corner block, transmissive, reflective, periodic, and user specified
                if (direction==ne) then
                    if (hydro_bound_type(2)==3) then
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,se,w_in)
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==se) then
                    if (hydro_bound_type(2)==3.and.hydro_bound_type(3)/=3) then
                        call extract_corner_hydro(blk_temp,sw,w_in)
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(2)/=3.and.hydro_bound_type(3)==3) then
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(2)==3.and.hydro_bound_type(3)==3) then
                        call extract_corner_hydro(blk_temp,nw,w_in)
                        w_out=w_in
                    end if
                else if (direction==sw) then
                    if (hydro_bound_type(3)==3) then
                        call extract_corner_hydro(blk_temp,ne,w_in)
                        w_out=w_in
                    else
                        call extract_corner_hydro(blk_temp,se,w_in)
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==nw) then
                    call extract_corner_hydro(blk_temp,se,w_in)
                    w_out=w_in
                end if
            case default
            end select
            if (direction==ne) then
                blk%w(blk_size_nx+1:blk_size_nx+2,blk_size_ny+1:blk_size_ny+2,1,1:5)=w_out
            else if (direction==se) then
                blk%w(blk_size_nx+1:blk_size_nx+2,-1:0,1,1:5)=w_out
            else if (direction==sw) then
                blk%w(-1:0,-1:0,1,1:5)=w_out
            else if (direction==nw) then
                blk%w(-1:0,blk_size_ny+1:blk_size_ny+2,1,1:5)=w_out
            end if
            nullify(blk_temp)
        end if
    end if
end subroutine communicate_hydro_corner_internal

subroutine assemble_communication_pattern()
    type(blockdef), pointer :: blk
    integer :: i,j,k
    if (nd==2) then
        processor%n_edge=0
        processor%n_corner=0
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            call count_external_edges(blk,processor%n_edge)
            call count_external_corners(blk,processor%n_corner)
            if (i/=np_nblk(rank+1)) then
                blk=>blk%blk_next
            end if
        end do
    end if
end subroutine assemble_communication_pattern

subroutine count_external_edges(blk,i)
    type(blockdef), pointer :: blk
    type(blockneighbour), pointer :: blk_nb
    integer :: i,j,neighbour_rank
    character(len=16) :: list(4)
    list=(/east,west,south,north/)
    do j=1,4
        nullify(blk_nb)
        call link_a_neighbour_2d(blk,list(j),blk_nb)
        if (associated(blk_nb)) then
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank/=rank) then
                i=i+1
            end if
        end if
    end do
end subroutine count_external_edges

subroutine count_external_corners(blk,i)
    type(blockdef), pointer :: blk
    type(blockneighbour), pointer :: blk_nb
    integer :: i,j,neighbour_rank
    character(len=16) :: list(4)
    list=(/se,sw,ne,nw/)
    do j=1,4
        nullify(blk_nb)
        call link_a_neighbour_2d(blk,list(j),blk_nb)
        if (associated(blk_nb)) then
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank/=rank) then
                i=i+1
            end if
        end if
    end do
end subroutine count_external_corners

subroutine communicate_hydro_external_2d()
    type(blockdef), pointer :: blk
    integer :: i,ierr,nreqs,ireq,msg_size,istart
    integer, allocatable :: reqs(:)
    real(8), dimension(:), allocatable :: msg
    nreqs=processor%n_edge+processor%n_corner
    msg_size=processor%n_edge*blk_xedge_size_hydro+processor%n_corner*blk_corner_size_hydro
    allocate(reqs(nreqs));reqs=-1;ireq=1
    allocate(msg(msg_size))
    istart=1
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        call communicate_hydro_edge_external_send(blk,east,reqs,ireq,istart,msg)
        call communicate_hydro_edge_external_send(blk,south,reqs,ireq,istart,msg)
        call communicate_hydro_edge_external_send(blk,west,reqs,ireq,istart,msg)
        call communicate_hydro_edge_external_send(blk,north,reqs,ireq,istart,msg)
        call communicate_hydro_corner_external_send(blk,se,reqs,ireq,istart,msg)
        call communicate_hydro_corner_external_send(blk,sw,reqs,ireq,istart,msg)
        call communicate_hydro_corner_external_send(blk,nw,reqs,ireq,istart,msg)
        call communicate_hydro_corner_external_send(blk,ne,reqs,ireq,istart,msg)
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        call communicate_hydro_edge_external_recv(blk,east)
        call communicate_hydro_edge_external_recv(blk,south)
        call communicate_hydro_edge_external_recv(blk,west)
        call communicate_hydro_edge_external_recv(blk,north)
        call communicate_hydro_corner_external_recv(blk,se)
        call communicate_hydro_corner_external_recv(blk,sw)
        call communicate_hydro_corner_external_recv(blk,nw)
        call communicate_hydro_corner_external_recv(blk,ne)
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    call mpi_waitall(nreqs,reqs,MPI_STATUSES_IGNORE,ierr)
    deallocate(reqs,msg)
end subroutine communicate_hydro_external_2d

subroutine communicate_hydro_edge_external_send(blk,direction,reqs,ireq,istart,msg)
    type(blockdef), pointer :: blk,blk_temp
    real(8), allocatable :: msg_out(:,:,:),msg(:)
    type(blockneighbour), pointer :: blk_nb
    integer, allocatable :: reqs(:)
    character(len=16) :: direction
    integer :: i,ireq,neighbour_rank,send_count,ierr,tag,stat(MPI_STATUS_SIZE),req,istart,iend
    call link_a_neighbour_2d(blk,direction,blk_nb)
    if (associated(blk_nb)) then
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            blk_temp=>blk_nb%blk
            if (direction==east) then
                send_count=2*blk_size_ny*5
                iend=istart+send_count-1
                allocate(msg_out(2,blk_size_ny,5))
                msg_out(1:2,1:blk_size_ny,1:5)=blk%w(blk_size_nx-1:blk_size_nx,1:blk_size_ny,1,1:5)
                msg(istart:iend)=reshape(msg_out,(/send_count/))
            else if (direction==south) then
                send_count=2*blk_size_nx*5
                iend=istart+send_count-1
                allocate(msg_out(blk_size_nx,2,5))
                msg_out(1:blk_size_nx,1:2,1:5)=blk%w(1:blk_size_nx,1:2,1,1:5)
                msg(istart:iend)=reshape(msg_out,(/send_count/))
            else if (direction==west) then
                send_count=2*blk_size_ny*5
                iend=istart+send_count-1
                allocate(msg_out(2,blk_size_ny,5))
                msg_out(1:2,1:blk_size_ny,1:5)=blk%w(1:2,1:blk_size_ny,1,1:5)
                msg(istart:iend)=reshape(msg_out,(/send_count/))
            else if (direction==north) then
                send_count=2*blk_size_nx*5
                iend=istart+send_count-1
                allocate(msg_out(blk_size_nx,2,5))
                msg_out(1:blk_size_nx,1:2,1:5)=blk%w(1:blk_size_nx,blk_size_ny-1:blk_size_ny,1,1:5)
                msg(istart:iend)=reshape(msg_out,(/send_count/))
            end if
            tag=hash_tag_send(blk,blk_temp,direction)
            call mpi_isend(msg(istart:iend),send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,reqs(ireq),ierr)
            ireq=ireq+1
            istart=iend+1
            deallocate(msg_out)
            nullify(blk_temp)
        end if
    end if
end subroutine communicate_hydro_edge_external_send

subroutine communicate_hydro_edge_external_recv(blk,direction)
    type(blockdef), pointer :: blk,blk_temp
    real(8), allocatable :: msg_in(:,:,:),msg(:)
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    integer :: i,j,neighbour_rank,recv_count,ierr,tag,stat(MPI_STATUS_SIZE),req
    call link_a_neighbour_2d(blk,direction,blk_nb)
    if (associated(blk_nb)) then
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            blk_temp=>blk_nb%blk
            if (direction==east.or.direction==west) then
                recv_count=2*blk_size_ny*5
                allocate(msg(recv_count),msg_in(2,blk_size_ny,5))
            else if (direction==south.or.direction==north) then
                recv_count=2*blk_size_nx*5
                allocate(msg(recv_count),msg_in(blk_size_nx,2,5))
            end if
            tag=hash_tag_recv(blk_temp,blk,direction)
            call mpi_recv(msg,recv_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
            msg_in=reshape(msg,shape(msg_in))
            if (direction==east) then
                blk%w(blk_size_nx+1:blk_size_nx+2,1:blk_size_ny,1,1:5)=msg_in(1:2,1:blk_size_ny,1:5)
            else if (direction==south) then
                blk%w(1:blk_size_nx,-1:0,1,1:5)=msg_in(1:blk_size_nx,1:2,1:5)
            else if (direction==west) then
                blk%w(-1:0,1:blk_size_ny,1,1:5)=msg_in(1:2,1:blk_size_ny,1:5)
            else if (direction==north) then
                blk%w(1:blk_size_nx,blk_size_ny+1:blk_size_ny+2,1,1:5)=msg_in(1:blk_size_nx,1:2,1:5)
            end if
            nullify(blk_temp)
            deallocate(msg_in,msg)
        end if
    end if
end subroutine communicate_hydro_edge_external_recv

subroutine communicate_hydro_corner_external_send(blk,direction,reqs,ireq,istart,msg)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    integer :: i,ireq,neighbour_rank,send_count,ierr,tag,stat(MPI_STATUS_SIZE),istart,iend
    integer, allocatable :: reqs(:)
    real(8), allocatable :: msg_out(:,:,:),msg(:)
    call link_a_neighbour_2d(blk,direction,blk_nb)
    if (associated(blk_nb)) then
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            blk_temp=>blk_nb%blk
            send_count=2*2*5
            iend=istart+send_count-1
            allocate(msg_out(2,2,5))
            if (direction==ne) then
                msg_out(1:2,1:2,1:5)=blk%w(blk_size_nx-1:blk_size_nx,blk_size_ny-1:blk_size_ny,1,1:5)
            else if (direction==se) then
                msg_out(1:2,1:2,1:5)=blk%w(blk_size_nx-1:blk_size_nx,1:2,1,1:5)
            else if (direction==sw) then
                msg_out(1:2,1:2,1:5)=blk%w(1:2,1:2,1,1:5)
            else if (direction==nw) then
                msg_out(1:2,1:2,1:5)=blk%w(1:2,blk_size_ny-1:blk_size_ny,1,1:5)
            end if
            tag=hash_tag_send(blk,blk_temp,direction)
            msg(istart:iend)=reshape(msg_out,(/send_count/))
            call mpi_isend(msg(istart:iend),send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,reqs(ireq),ierr)
            istart=iend+1
            ireq=ireq+1
            nullify(blk_temp)
            deallocate(msg_out)
        end if
    end if
end subroutine communicate_hydro_corner_external_send

subroutine communicate_hydro_corner_external_recv(blk,direction)
    !all cases of intercore corner communications are in this subroutine
    !direction indicates the direction of the neighbour block
    !for boundary blocks, the guard corner may be mirrored
    !reflect means reflective boundary
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    integer :: i,j,neighbour_rank,recv_count,ierr,tag,req,stat(MPI_STATUS_SIZE)
    real(8), allocatable :: msg_in(:,:,:),msg(:)
    real(8), dimension(2,2,5) :: w_in,u_in,w_out,u_out
    real(8), dimension(2,2) :: temp_in,egv_in,temp_out,egv_out
    character(len=1) :: mirror
    logical :: reflect
    call link_a_neighbour_2d(blk,direction,blk_nb)
    if (associated(blk_nb)) then
        neighbour_rank=blk_nb%cpu_rank
        if (neighbour_rank/=rank) then
            blk_temp=>blk_nb%blk
            recv_count=2*2*5
            allocate(msg_in(2,2,5),msg(recv_count))
            tag=hash_tag_recv(blk_temp,blk,direction)
            call mpi_recv(msg,recv_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
            msg_in=reshape(msg,shape(msg_in))
            w_in=msg_in(1:2,1:2,1:5)
            select case (blk%loc_type)
            case (1)
                if (direction==ne) then
                    if (hydro_bound_type(4)==3) then
                        w_out=w_in
                    else
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==se) then
                    w_out=w_in
                else if (direction==sw) then
                    if (hydro_bound_type(1)==3) then
                        w_out=w_in
                    else
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==nw) then
                    if (hydro_bound_type(1)==3.and.hydro_bound_type(4)/=3) then
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(1)/=3.and.hydro_bound_type(4)==3) then
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(1)==3.and.hydro_bound_type(4)==3) then
                        w_out=w_in
                    end if
                end if
            case (2)
                if (direction==ne.or.direction==nw) then
                    select case (hydro_bound_type(4))
                    case (1)
                        mirror='y';reflect=.false.
                        call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                    case (2)
                        mirror='y';reflect=.true.
                        call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                    case (3)
                        w_out=w_in
                    end select
                else
                    w_out=w_in
                end if
            case (3)
                if (direction==ne) then
                    if (hydro_bound_type(2)==3.and.hydro_bound_type(4)/=3) then
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(2)/=3.and.hydro_bound_type(4)==3) then
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(2)==3.and.hydro_bound_type(4)==3) then
                        w_out=w_in
                    end if
                else if (direction==se) then
                    if (hydro_bound_type(2)==3) then
                        w_out=w_in
                    else
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==sw) then
                    w_out=w_in
                else if (direction==nw) then
                    if (hydro_bound_type(4)==3) then
                        w_out=w_in
                    else
                        if (hydro_bound_type(4)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(4)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                end if
            case (4)
                if (direction==nw.or.direction==sw) then
                    select case (hydro_bound_type(1))
                    case (1)
                        mirror='x';reflect=.false.
                        call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                    case (2)
                        mirror='x';reflect=.true.
                        call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                    case (3)
                        w_out=w_in
                    end select
                else
                    w_out=w_in
                end if
            case (5)
                !inside, no mirror or transform
                w_out=w_in
            case (6)
                if (direction==ne.or.direction==se) then
                    select case (hydro_bound_type(2))
                    case (1)
                        mirror='x';reflect=.false.
                        call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                    case (2)
                        mirror='x';reflect=.false.
                        call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                    case (3)
                        w_out=w_in
                    end select
                else
                    w_out=w_in
                end if
            case (7)
                if (direction==ne) then
                    w_out=w_in
                else if (direction==se) then
                    if (hydro_bound_type(3)==3) then
                        w_out=w_in
                    else
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==sw) then
                    if (hydro_bound_type(1)==3.and.hydro_bound_type(3)/=3) then
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(1)/=3.and.hydro_bound_type(3)==3) then
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(1)==3.and.hydro_bound_type(3)==3) then
                        w_out=w_in
                    end if
                else if (direction==nw) then
                    if (hydro_bound_type(1)==3) then
                        w_out=w_in
                    else
                        if (hydro_bound_type(1)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(1)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                end if
            case (8)
                if (direction==sw.or.direction==se) then
                    select case (hydro_bound_type(3))
                    case (1)
                        mirror='y';reflect=.false.
                        call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                    case (2)
                        mirror='y';reflect=.true.
                        call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                    case (3)
                        w_out=w_in
                    end select
                else
                    w_out=w_in
                end if
            case (9)
                if (direction==ne) then
                    if (hydro_bound_type(2)==3) then
                        w_out=w_in
                    else
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==se) then
                    if (hydro_bound_type(2)==3.and.hydro_bound_type(3)/=3) then
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(2)/=3.and.hydro_bound_type(3)==3) then
                        if (hydro_bound_type(2)==1) then
                            mirror='x';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(2)==2) then
                            mirror='x';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    else if (hydro_bound_type(2)==3.and.hydro_bound_type(3)==3) then
                        w_out=w_in
                    end if
                else if (direction==sw) then
                    if (hydro_bound_type(3)==3) then
                        w_out=w_in
                    else
                        if (hydro_bound_type(3)==1) then
                            mirror='y';reflect=.false.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        else if (hydro_bound_type(3)==2) then
                            mirror='y';reflect=.true.
                            call corner_mirror_reflect_hydro(w_in,mirror,reflect,w_out)
                        end if
                    end if
                else if (direction==nw) then
                    w_out=w_in
                end if
            case default
            end select
            if (direction==ne) then
                blk%w(blk_size_nx+1:blk_size_nx+2,blk_size_ny+1:blk_size_ny+2,1,1:5)=w_out
            else if (direction==se) then
                blk%w(blk_size_nx+1:blk_size_nx+2,-1:0,1,1:5)=w_out
            else if (direction==sw) then
                blk%w(-1:0,-1:0,1,1:5)=w_out
            else if (direction==nw) then
                blk%w(-1:0,blk_size_ny+1:blk_size_ny+2,1,1:5)=w_out
            end if
            nullify(blk_temp,blk_nb)
            deallocate(msg_in)
        end if
    end if
end subroutine communicate_hydro_corner_external_recv

subroutine communicate_blocks_fld_1d()
    call communicate_blocks_fld_internal_1d()
    call communicate_blocks_fld_external_1d()
end subroutine communicate_blocks_fld_1d

subroutine communicate_blocks_fld_internal_1d()
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    integer :: i,neighbour_rank
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        if (associated(blk%nb_l)) then
            blk_nb=>blk%nb_l
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank==rank) then
                blk_temp=>blk_nb%blk
                blk%Erad_pre=blk_temp%Erad(blk_size_nx,1,1)
                blk%Erad_pre_int=blk_temp%Erad_int(blk_size_nx,1,1)
                blk%sigma_rosseland_pre=blk_temp%sigma_rosseland(blk_size_nx,1,1)
                blk%nb_coor%xl(1)=blk_temp%x_center(blk_size_nx)
            end if
        end if
        if (associated(blk%nb_r)) then
            blk_nb=>blk%nb_r
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank==rank) then
                blk_temp=>blk_nb%blk
                blk%Erad_next=blk_temp%Erad(1,1,1)
                blk%Erad_next_int=blk_temp%Erad_int(1,1,1)
                blk%sigma_rosseland_next=blk_temp%sigma_rosseland(1,1,1)
                blk%nb_coor%xu(1)=blk_temp%x_center(1)
            end if
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    nullify(blk)
end subroutine communicate_blocks_fld_internal_1d

subroutine communicate_blocks_fld_external_1d()
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    real(8) :: msg_out(4),msg_in(4)
    integer :: i,neighbour_rank,req,ierr,stat(MPI_STATUS_SIZE)
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        if (associated(blk%nb_l)) then
            blk_nb=>blk%nb_l
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank/=rank) then
                msg_out(1)=blk%Erad(1,1,1)
                msg_out(2)=blk%Erad_int(1,1,1)
                msg_out(3)=blk%sigma_rosseland(1,1,1)
                msg_out(4)=blk%x_center(1)
                call mpi_isend(msg_out,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,req,ierr)
                call mpi_recv(msg_in,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,stat,ierr)
                blk%Erad_pre=msg_in(1)
                blk%Erad_pre_int=msg_in(2)
                blk%sigma_rosseland_pre=msg_in(3)
                blk%nb_coor%xl(1)=msg_in(4)
            end if
        end if
        if (associated(blk%nb_r)) then
            blk_nb=>blk%nb_r
            neighbour_rank=blk_nb%cpu_rank
            if (neighbour_rank/=rank) then
                msg_out(1)=blk%Erad(blk_size_nx,1,1)
                msg_out(2)=blk%Erad_int(blk_size_nx,1,1)
                msg_out(3)=blk%sigma_rosseland(blk_size_nx,1,1)
                msg_out(4)=blk%x_center(blk_size_nx)
                call mpi_isend(msg_out,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,req,ierr)
                call mpi_recv(msg_in,4,MPI_REAL8,neighbour_rank,1,MPI_COMM_WORLD,stat,ierr)
                blk%Erad_next=msg_in(1)
                blk%Erad_next_int=msg_in(2)
                blk%sigma_rosseland_next=msg_in(3)
                blk%nb_coor%xu(1)=msg_in(4)
            end if
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    nullify(blk)
end subroutine communicate_blocks_fld_external_1d

subroutine communicate_blocks_2d_fld()
    call communicate_east_west_fld()
    call communicate_south_north_fld()
end subroutine communicate_blocks_2d_fld

subroutine communicate_east_fld(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    integer :: neighbour_rank,send_count,recv_count,ierr,tag,req,stat(MPI_STATUS_SIZE)
    real(8), allocatable :: msg_in(:,:),msg_out(:,:)
    blk_nb=>blk%nb_e
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank==rank) then
        blk_temp=>blk_nb%blk
        blk%Erad_xu=blk_temp%Erad(1,1:blk_size_ny,1)
        blk%sigma_rosseland_xu=blk_temp%sigma_rosseland(1,1:blk_size_ny,1)
        blk%dx_xu=blk_temp%dxyz(1)
    else
        send_count=2*blk_size_ny
        recv_count=send_count
        allocate(msg_in(blk_size_ny,2),msg_out(blk_size_ny,2))
        msg_out(:,1)=blk%Erad(blk_size_nx,1:blk_size_ny,1)
        msg_out(:,2)=blk%sigma_rosseland(blk_size_nx,1:blk_size_ny,1)
        tag=blk%blk_id*10;call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        tag=blk%blk_id*10+1;call mpi_isend(blk%dxyz(1),1,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        tag=blk_nb%blk%blk_id*10;call mpi_recv(msg_in,recv_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        tag=blk_nb%blk%blk_id*10+1;call mpi_recv(blk%dx_xu,1,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        blk%Erad_xu=msg_in(:,1)
        blk%sigma_rosseland_xu=msg_in(:,2)
        deallocate(msg_in,msg_out)
    end if
    nullify(blk_nb,blk_temp)
end subroutine communicate_east_fld

subroutine communicate_west_fld(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    integer :: neighbour_rank,send_count,recv_count,ierr,tag,req,stat(MPI_STATUS_SIZE)
    real(8), allocatable :: msg_in(:,:),msg_out(:,:)
    blk_nb=>blk%nb_w
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank==rank) then
        blk_temp=>blk_nb%blk
        blk%Erad_xl=blk_temp%Erad(blk_size_nx,1:blk_size_ny,1)
        blk%sigma_rosseland_xl=blk_temp%sigma_rosseland(blk_size_nx,1:blk_size_ny,1)
        blk%dx_xl=blk_temp%dxyz(1)
    else
        send_count=2*blk_size_ny
        recv_count=send_count
        allocate(msg_in(blk_size_ny,2),msg_out(blk_size_ny,2))
        tag=blk_nb%blk%blk_id*10;call mpi_recv(msg_in,recv_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        tag=blk_nb%blk%blk_id*10+1;call mpi_recv(blk%dx_xl,1,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        msg_out(:,1)=blk%Erad(1,1:blk_size_ny,1)
        msg_out(:,2)=blk%sigma_rosseland(1,1:blk_size_ny,1)
        tag=blk%blk_id*10;call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        tag=blk%blk_id*10+1;call mpi_isend(blk%dxyz(1),1,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        blk%Erad_xl=msg_in(:,1)
        blk%sigma_rosseland_xl=msg_in(:,2)
        deallocate(msg_in,msg_out)
    end if
    nullify(blk_nb,blk_temp)
end subroutine communicate_west_fld

subroutine communicate_east_west_fld()
    type(blockdef), pointer :: blk_traversal
    integer :: i,ierr
    blk_traversal=>blk_processor_head
    do i=1,np_nblk(rank+1)
        select case (blk_traversal%loc_type)
        case (1)
            call communicate_east_fld(blk_traversal)
        case (2)
            call communicate_east_fld(blk_traversal)
            call communicate_west_fld(blk_traversal)
        case (3)
            call communicate_west_fld(blk_traversal)
        case (4)
            call communicate_east_fld(blk_traversal)
        case (5)
            call communicate_east_fld(blk_traversal)
            call communicate_west_fld(blk_traversal)
        case (6)
            call communicate_west_fld(blk_traversal)
        case (7)
            call communicate_east_fld(blk_traversal)
        case (8)
            call communicate_east_fld(blk_traversal)
            call communicate_west_fld(blk_traversal)
        case (9)
            call communicate_west_fld(blk_traversal)
        end select
        if (i/=np_nblk(rank+1)) then
            blk_traversal=>blk_traversal%blk_next
        end if
    end do
    nullify(blk_traversal)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
end subroutine communicate_east_west_fld

subroutine communicate_south_fld(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    integer :: neighbour_rank,send_count,recv_count,ierr,tag,req,stat(MPI_STATUS_SIZE)
    real(8), allocatable :: msg_in(:,:),msg_out(:,:)
    blk_nb=>blk%nb_s
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank==rank) then
        blk_temp=>blk_nb%blk
        blk%Erad_yl=blk_temp%Erad(1:blk_size_nx,blk_size_ny,1)
        blk%sigma_rosseland_yl=blk_temp%sigma_rosseland(1:blk_size_nx,blk_size_ny,1)
        blk%dx_yl=blk_temp%dxyz(1)
    else
        send_count=2*blk_size_nx
        recv_count=send_count
        allocate(msg_in(blk_size_nx,2),msg_out(blk_size_nx,2))
        msg_out(1:blk_size_nx,1)=blk%Erad(1:blk_size_nx,1,1)
        msg_out(1:blk_size_nx,2)=blk%sigma_rosseland(1:blk_size_nx,1,1)
        tag=blk%blk_id*10;call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        tag=blk%blk_id*10+1;call mpi_isend(blk%dxyz(1),1,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        tag=blk_nb%blk%blk_id*10;call mpi_recv(msg_in,recv_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        tag=blk_nb%blk%blk_id*10+1;call mpi_recv(blk%dx_yl,1,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        blk%Erad_yl=msg_in(1:blk_size_nx,1)
        blk%sigma_rosseland_yl=msg_in(1:blk_size_nx,2)
        call mpi_wait(req,stat,ierr)
        deallocate(msg_in,msg_out)
    end if
    nullify(blk_nb,blk_temp)
end subroutine communicate_south_fld

subroutine communicate_north_fld(blk)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    integer :: neighbour_rank,send_count,recv_count,ierr,tag,req,stat(MPI_STATUS_SIZE)
    real(8), allocatable :: msg_in(:,:),msg_out(:,:)
    blk_nb=>blk%nb_n
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank==rank) then
        blk_temp=>blk_nb%blk
        blk%Erad_yu=blk_temp%Erad(1:blk_size_nx,1,1)
        blk%sigma_rosseland_yu=blk_temp%sigma_rosseland(1:blk_size_nx,1,1)
        blk%dx_yu=blk_temp%dxyz(1)
    else
        send_count=2*blk_size_nx
        recv_count=send_count
        allocate(msg_in(blk_size_nx,2),msg_out(blk_size_nx,2))
        tag=blk_nb%blk%blk_id*10;call mpi_recv(msg_in,recv_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        tag=blk_nb%blk%blk_id*10+1;call mpi_recv(blk%dx_yu,1,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        msg_out(1:blk_size_nx,1)=blk%Erad(1:blk_size_nx,blk_size_ny,1)
        msg_out(1:blk_size_nx,2)=blk%sigma_rosseland(1:blk_size_nx,blk_size_ny,1)
        tag=blk%blk_id*10;call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        tag=blk%blk_id*10+1;call mpi_isend(blk%dxyz(1),1,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        blk%Erad_yu=msg_in(1:blk_size_nx,1)
        blk%sigma_rosseland_yu=msg_in(1:blk_size_nx,2)
        call mpi_wait(req,stat,ierr)
        deallocate(msg_in,msg_out)
    end if
    nullify(blk_nb,blk_temp)
end subroutine communicate_north_fld

subroutine communicate_south_north_fld()
    type(blockdef), pointer :: blk_traversal
    integer :: i,ierr
    blk_traversal=>blk_processor_head
    do i=1,np_nblk(rank+1)
        select case (blk_traversal%loc_type)
        case (1)
            call communicate_south_fld(blk_traversal)
        case (2)
            call communicate_south_fld(blk_traversal)
        case (3)
            call communicate_south_fld(blk_traversal)
        case (4)
            call communicate_south_fld(blk_traversal)
            call communicate_north_fld(blk_traversal)
        case (5)
            call communicate_south_fld(blk_traversal)
            call communicate_north_fld(blk_traversal)
        case (6)
            call communicate_south_fld(blk_traversal)
            call communicate_north_fld(blk_traversal)
        case (7)
            call communicate_north_fld(blk_traversal)
        case (8)
            call communicate_north_fld(blk_traversal)
        case (9)
            call communicate_north_fld(blk_traversal)
        end select
        if (i/=np_nblk(rank+1)) then
            blk_traversal=>blk_traversal%blk_next
        end if
    end do
    nullify(blk_traversal)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
end subroutine communicate_south_north_fld

subroutine block_restriction(blk)
    type(blockdef), pointer :: blk,blk_xl,blk_xu
    real(8) :: u(5),w(5),temp,egv,Erad
    integer :: i,j,k
    if (nd==1) then
        blk%on_processor=.true.
        blk_xl=>blk%blk_xl
        blk_xu=>blk%blk_xu
        do i=blk_xlb+1,blk_size_nx/2
            u=(blk_xl%u(2*i-1,1,1,1:5)+blk_xl%u(2*i,1,1,1:5))/2
#if         ieos==1
            call utow(u,w)
            temp=calculate_temp_from_p_rho(w(5),w(1))
            egv=w(5)/(gamma_gas-1)
#elif       ieos==2
            call eos_hllc_analytic_utow(u,w,temp,egv)
#endif
            blk%u(i,1,1,1:5)=u
            blk%w(i,1,1,1:5)=w
            blk%temp(i,1,1)=temp
            blk%egv(i,1,1)=egv
            if (iradiation/=0) then
                blk%Erad(i,1,1)=(blk_xl%Erad(2*i-1,1,1)+blk_xl%Erad(2*i,1,1))/2
            end if
        end do
        do i=blk_size_nx/2+1,blk_xub-1
            u=(blk_xu%u(2*(i-blk_size_nx/2)-1,1,1,1:5)+blk_xu%u(2*(i-blk_size_nx/2),1,1,1:5))/2
#if         ieos==1
            call utow(u,w)
            temp=calculate_temp_from_p_rho(w(5),w(1))
            egv=w(5)/(gamma_gas-1)
#elif       ieos==2
            call eos_hllc_analytic_utow(u,w,temp,egv)
#endif
            blk%u(i,1,1,1:5)=u
            blk%w(i,1,1,1:5)=w
            blk%temp(i,1,1)=temp
            blk%egv(i,1,1)=egv
            if (iradiation/=0) then
                blk%Erad(i,1,1)=(blk_xu%Erad(2*(i-blk_size_nx/2)-1,1,1)+blk_xu%Erad(2*(i-blk_size_nx/2),1,1))/2
            end if
        end do
    else if (nd==2) then
    end if
end subroutine block_restriction

subroutine block_prolongation(blk)
    !blk is a parent block, prolongate it and store the data to its children
    type(blockdef), pointer :: blk
    real(8) :: diff_l,diff_r,diff_c,diff,u(2,5),w(2,5),temp(2),egv(2),u1(5),w1(5)
    integer :: i,j,k
    if (nd==1) then
        blk%blk_xl%on_processor=.true.
        blk%blk_xu%on_processor=.true.
        do i=blk_xlb+1,blk_size_nx/2+1
            do j=1,5
                diff_l=blk%u(i,1,1,j)-blk%u(i-1,1,1,j)
                diff_r=blk%u(i+1,1,1,j)-blk%u(i,1,1,j)
                diff_c=(blk%u(i+1,1,1,j)-blk%u(i-1,1,1,j))/2
                if (diff_c/=zero) then
                    diff=abs(diff_c)/diff_c*min(abs(diff_l),abs(diff_r),abs(diff_c))
                else
                    diff=0d0
                end if
                diff=0d0
                u(1,j)=blk%u(i,1,1,j)-diff/4
                u(2,j)=blk%u(i,1,1,j)+diff/4
            end do
#if         ieos==1
            do j=1,2
                u1=u(j,1:5)
                call utow(u1,w1)
                w(j,1:5)=w1
                temp(j)=calculate_temp_from_p_rho(w1(5),w1(1))
                egv(j)=w1(5)/(gamma_gas-one)
            end do
#elif       ieos==2
            do j=1,2
                u1=u(j,1:5)
                call eos_hllc_analytic_utow(u1,w1,temp(j),egv(j))
                w(j,1:5)=w1
            end do
#endif
            blk%blk_xl%u(2*i-1,1,1,1:5)=u(1,1:5)
            blk%blk_xl%w(2*i-1,1,1,1:5)=w(1,1:5)
            blk%blk_xl%temp(2*i-1,1,1)=temp(1)
            blk%blk_xl%egv(2*i-1,1,1)=egv(1)
            blk%blk_xl%u(2*i,1,1,1:5)=u(2,1:5)
            blk%blk_xl%w(2*i,1,1,1:5)=w(2,1:5)
            blk%blk_xl%temp(2*i,1,1)=temp(2)
            blk%blk_xl%egv(2*i,1,1)=egv(2)
            if (iradiation/=0) then
                blk%blk_xl%Erad(2*i-1,1,1)=blk%Erad(i,1,1)
                blk%blk_xl%Erad(2*i,1,1)=blk%Erad(i,1,1)
            end if
        end do
        do i=blk_size_nx/2,blk_xub-1
            do j=1,5
                diff_l=blk%u(i,1,1,j)-blk%u(i-1,1,1,j)
                diff_r=blk%u(i+1,1,1,j)-blk%u(i,1,1,j)
                diff_c=(blk%u(i+1,1,1,j)-blk%u(i-1,1,1,j))/2
                if (diff_c/=zero) then
                    diff=abs(diff_c)/diff_c*min(abs(diff_l),abs(diff_r),abs(diff_c))
                else
                    diff=0d0
                end if
                diff=0d0
                u(1,j)=blk%u(i,1,1,j)-diff/4
                u(2,j)=blk%u(i,1,1,j)+diff/4
            end do
#if         ieos==1
            do j=1,2
                u1=u(j,1:5)
                call utow(u1,w1)
                w(j,1:5)=w1
                temp(j)=calculate_temp_from_p_rho(w1(5),w1(1))
                egv(j)=w1(5)/(gamma_gas-one)
            end do
#elif       ieos==2
            do j=1,2
                u1=u(j,1:5)
                call eos_hllc_analytic_utow(u1,w1,temp(j),egv(j))
                w(j,1:5)=w1
            end do
#endif
            blk%blk_xu%u(2*(i-blk_size_nx/2)-1,1,1,1:5)=u(1,1:5)
            blk%blk_xu%w(2*(i-blk_size_nx/2)-1,1,1,1:5)=w(1,1:5)
            blk%blk_xu%temp(2*(i-blk_size_nx/2)-1,1,1)=temp(1)
            blk%blk_xu%egv(2*(i-blk_size_nx/2)-1,1,1)=egv(1)
            blk%blk_xu%u(2*(i-blk_size_nx/2),1,1,1:5)=u(2,1:5)
            blk%blk_xu%w(2*(i-blk_size_nx/2),1,1,1:5)=w(2,1:5)
            blk%blk_xu%temp(2*(i-blk_size_nx/2),1,1)=temp(2)
            blk%blk_xu%egv(2*(i-blk_size_nx/2),1,1)=egv(2)
            if (iradiation/=0) then
                blk%blk_xu%Erad(2*(i-blk_size_nx/2)-1,1,1)=blk%Erad(i,1,1)
                blk%blk_xu%Erad(2*(i-blk_size_nx/2),1,1)=blk%Erad(i,1,1)
            end if
        end do
    else if (nd==2) then
    end if
end subroutine block_prolongation

subroutine reset_smr_processor_variables()
    !upon restart, the smr structure may change
    processor%smr_grow=.false.
    processor%smr_derefine=.false.
    processor%smr_ngrow=0
    processor%smr_nderefine=0
    if (.not.allocated(processor%mpi_logical)) allocate(processor%mpi_logical(np))
    processor%mpi_logical=.false.
    if (.not.allocated(processor%mpi_processor_oper)) allocate(processor%mpi_processor_oper(np))
    processor%mpi_processor_oper=0
    if (allocated(processor%grow_keys)) deallocate(processor%grow_keys)
    if (allocated(processor%derefine_keys)) deallocate(processor%derefine_keys)
    if (allocated(processor%mpi_world_keys)) deallocate(processor%mpi_world_keys)
end subroutine reset_smr_processor_variables

subroutine reset_amr_processor_variables()
    processor%amr_grow=.false.
    processor%amr_derefine=.false.
    processor%amr_ngrow=0
    processor%amr_nderefine=0
    if (.not.allocated(processor%mpi_logical)) allocate(processor%mpi_logical(np))
    processor%mpi_logical=.false.
    if (.not.allocated(processor%mpi_processor_oper)) allocate(processor%mpi_processor_oper(np))
    processor%mpi_processor_oper=0
    if (allocated(processor%grow_keys)) deallocate(processor%grow_keys)
    if (allocated(processor%derefine_keys)) deallocate(processor%derefine_keys)
    if (allocated(processor%mpi_world_keys)) deallocate(processor%mpi_world_keys)
end subroutine reset_amr_processor_variables

subroutine record_processor_smr_grow_order(blk)
    !record the keys of the blocks that are refined
    type(blockdef), pointer :: blk
    integer, dimension(:,:), allocatable :: temp_keys
end subroutine record_processor_smr_grow_order

subroutine record_processor_amr_grow_order(blk)
    !record the keys of the blocks that are refined
    type(blockdef), pointer :: blk
    integer, dimension(:,:), allocatable :: temp_keys
    processor%amr_grow=.true.
    processor%amr_ngrow=processor%amr_ngrow+1
    if (processor%amr_ngrow==1) then
        !the first grown key
        allocate(processor%grow_keys(3,1))
        processor%grow_keys(1:3,1)=blk%key
    else
        !add another key to the existing list
        allocate(temp_keys(3,processor%amr_ngrow))
        temp_keys(1:3,1:processor%amr_ngrow-1)=processor%grow_keys
        temp_keys(1:3,processor%amr_ngrow)=blk%key
        deallocate(processor%grow_keys)
        allocate(processor%grow_keys(3,processor%amr_ngrow))
        processor%grow_keys=temp_keys
        deallocate(temp_keys)
    end if
end subroutine record_processor_amr_grow_order

subroutine record_processor_amr_derefine_order(blk)
    !record the keys of the blocks that are derefined
    type(blockdef), pointer :: blk
    integer, dimension(:,:), allocatable :: temp_keys
    processor%amr_derefine=.true.
    processor%amr_nderefine=processor%amr_nderefine+1
    if (processor%amr_nderefine==1) then
        allocate(processor%derefine_keys(3,1))
        processor%derefine_keys(1:3,1)=blk%key
    else
        allocate(temp_keys(3,processor%amr_nderefine))
        temp_keys(1:3,1:processor%amr_nderefine-1)=processor%derefine_keys
        temp_keys(1:3,processor%amr_nderefine)=blk%key
        deallocate(processor%derefine_keys)
        allocate(processor%derefine_keys(3,processor%amr_nderefine))
        processor%derefine_keys=temp_keys
        deallocate(temp_keys)
    end if
end subroutine record_processor_amr_derefine_order

subroutine sync_the_tree_grow()
    type(blockdef), pointer :: blk
    integer :: ierr,i,j,n_oper,reqs,p_start,send_count,receive_count,key(3)
    integer, dimension(:,:), allocatable :: keys
    call mpi_gather(processor%amr_grow,1,MPI_LOGICAL,processor%mpi_logical,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(processor%mpi_logical,np,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_gather(processor%amr_ngrow,1,MPI_INTEGER,processor%mpi_processor_oper,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(processor%mpi_processor_oper,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    np_nblk=np_nblk+processor%mpi_processor_oper
    n_oper=sum(processor%mpi_processor_oper)
    if (n_oper>0) then
        reqs=0
        allocate(processor%mpi_world_keys(3,n_oper))
        processor%mpi_world_keys=0
        if (processor%amr_ngrow>0) then
            send_count=3*processor%amr_ngrow
            call mpi_isend(processor%grow_keys,send_count,MPI_INTEGER,0,1,MPI_COMM_WORLD,reqs,ierr)
        end if
        if (rank==0) then
            !stack all the keys of blocks that will be refined
            do i=0,np-1
                if (i>0) then
                    p_start=sum(processor%mpi_processor_oper(1:i))+1
                else
                    p_start=1
                end if
                receive_count=3*processor%mpi_processor_oper(i+1)
                if (receive_count>0) then
                    call mpi_irecv(processor%mpi_world_keys(1,p_start),receive_count,MPI_INTEGER,i,1,MPI_COMM_WORLD,reqs,ierr)
                end if
            end do
        end if
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        send_count=size(processor%mpi_world_keys)
        call mpi_bcast(processor%mpi_world_keys,send_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        do i=1,n_oper
            key=processor%mpi_world_keys(1:3,i)
            call block_key_to_pointer(key,blk)
            call grow_tree_node_1d_conditional(blk)
        end do
    end if
end subroutine sync_the_tree_grow

subroutine label_discontinuities()
end subroutine label_discontinuities

subroutine sync_the_tree_derefine()
    type(blockdef), pointer :: blk
    integer :: ierr,i,j,n_oper,reqs,p_start,send_count,receive_count,key(3)
    integer, dimension(:,:), allocatable :: keys
    call mpi_gather(processor%amr_derefine,1,MPI_LOGICAL,processor%mpi_logical,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(processor%mpi_logical,np,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_gather(processor%amr_nderefine,1,MPI_INTEGER,processor%mpi_processor_oper,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(processor%mpi_processor_oper,np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    np_nblk=np_nblk-processor%mpi_processor_oper
    n_oper=sum(processor%mpi_processor_oper)
    if (n_oper>0) then
        reqs=0
        allocate(processor%mpi_world_keys(3,n_oper))
        processor%mpi_world_keys=0
        if (processor%amr_nderefine>0) then
            send_count=3*processor%amr_nderefine
            call mpi_isend(processor%derefine_keys,send_count,MPI_INTEGER,0,1,MPI_COMM_WORLD,reqs,ierr)
        end if
        if (rank==0) then
            !gather all the keys of blocks that will be derefined
            do i=0,np-1
                if (i>0) then
                    p_start=sum(processor%mpi_processor_oper(1:i))+1
                else
                    p_start=1
                end if
                receive_count=3*processor%mpi_processor_oper(i+1)
                if (receive_count>0) then
                    call mpi_irecv(processor%mpi_world_keys(1,p_start),receive_count,MPI_INTEGER,i,1,MPI_COMM_WORLD,reqs,ierr)
                end if
            end do
        end if
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        send_count=size(processor%mpi_world_keys)
        call mpi_bcast(processor%mpi_world_keys,send_count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        do i=1,n_oper
            key=processor%mpi_world_keys(1:3,i)
            call block_key_to_pointer(key,blk)
            call derefine_amr_node_conditional(blk)
        end do
    end if
end subroutine sync_the_tree_derefine

subroutine sync_renumber_reset()
    call sync_the_tree_grow()
    call renumber_domain_blocks()
    call reset_amr_processor_variables()
end subroutine sync_renumber_reset

subroutine grow_amr_tree()
    type(blockdef), pointer :: blk
    integer :: i,j,ierr
    if (refine_type=='adaptive'.or.refine_type=='mixed') then
        if (nd==1) then
            call reset_amr_processor_variables()
            blk=>blk_processor_head
            do i=1,np_nblk(rank+1)
                call grow_amr_node(blk)
                blk=>blk%blk_next
            end do
            call sync_renumber_reset()
            call label_discontinuities()
            do i=1,max_refine_level-1
                blk=>blk_processor_head
                do while (.not.associated(blk,blk_processor_tail))
                    call add_amr_buffer_blocks(blk)
                    blk=>blk%blk_next
                end do
                call add_amr_buffer_blocks(blk)
            end do
            call sync_renumber_reset()
            call examine_derefine()
            blk=>blk_processor_head
            do while (.not.associated(blk,blk_processor_tail))
                if (blk%level>0) then
                    call derefine_amr_node(blk)
                end if
                if (associated(blk,blk_processor_tail)) then
                    exit
                else
                    blk=>blk%blk_next
                end if
            end do
            nullify(blk)
            call sync_the_tree_derefine()
            call amr_load_balancer()
            call identify_block_spatial_type()
            call link_neighbour_blocks()
            call mpi_barrier(MPI_COMM_WORLD,ierr)
            call communicate_hydro()
        else if (nd==2) then
        end if
    end if
end subroutine grow_amr_tree

subroutine amr_load_balancer()
    type(blockdef), pointer :: blk
    integer :: blks_per_core,remainder,i,j,blk_id,former_rank,new_rank,ierr
    integer, dimension(:), allocatable :: np_nblk_old
    integer, dimension(:,:), allocatable :: stats
    if (nd==1) then
        call renumber_domain_blocks()
        allocate(np_nblk_old(np))
        np_nblk_old=np_nblk
        blks_per_core=nblk_total/np
        np_nblk=blks_per_core
        remainder=mod(nblk_total,np)
        do i=1,remainder
            np_nblk(i)=np_nblk(i)+1
        end do
        !send the block data to the another processor if neccessary
        blk=>blk_processor_head
        do i=1,np_nblk_old(rank+1)
            blk_id=blk%blk_id
            call block_id_to_processor_rank(np_nblk,blk_id,new_rank)
            if (rank/=new_rank) then
                call send_block(blk,new_rank)
            end if
            blk=>blk%blk_next
        end do
        !locate the new blk_processor_head and blk_processor_tail
        blk=>blk_head
        do i=1,rank+1
            blk_processor_head=>blk
            do j=1,np_nblk(i)-1
                blk=>blk%blk_next
            end do
            blk_processor_tail=>blk
            if (associated(blk%blk_next)) then
                blk=>blk%blk_next
            end if
        end do
        !receive the block data from another processor
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            blk_id=blk%blk_id
            call block_id_to_processor_rank(np_nblk_old,blk_id,former_rank)
            if (rank/=former_rank) then
                if (allocated(blk%w)) then
                    call receive_block(blk,former_rank)
                else
                    call allocate_block_heavy_data(blk)
                    call receive_block(blk,former_rank)
                end if
            end if
            blk=>blk%blk_next
        end do
    else if (nd==2) then
    end if
    blk=>blk_head
    do i=1,nblk_total
        blk%on_processor=.false.
        blk=>blk%blk_next
    end do
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        blk%on_processor=.true.
        blk=>blk%blk_next
    end do
    deallocate(np_nblk_old)
end subroutine amr_load_balancer

subroutine restart_smr_load_balancer()
    type(blockdef), pointer :: blk
    integer :: blks_per_core,remainder,i,j,blk_id,former_rank,new_rank,ierr
    integer, dimension(:), allocatable :: np_nblk_old
    integer, dimension(:,:), allocatable :: stats
    if (nd==1) then
        call renumber_domain_blocks()
        allocate(np_nblk_old(np))
        np_nblk_old=np_nblk
        blks_per_core=nblk_total/np
        np_nblk=blks_per_core
        remainder=mod(nblk_total,np)
        do i=1,remainder
            np_nblk(i)=np_nblk(i)+1
        end do
        !send the block data to the another processor if neccessary
        blk=>blk_processor_head
        do i=1,np_nblk_old(rank+1)
            blk_id=blk%blk_id
            call block_id_to_processor_rank(np_nblk,blk_id,new_rank)
            if (rank/=new_rank) then
                call send_block(blk,new_rank)
            end if
            blk=>blk%blk_next
        end do
        !locate the new blk_processor_head and blk_processor_tail
        blk=>blk_head
        do i=1,rank+1
            blk_processor_head=>blk
            do j=1,np_nblk(i)-1
                blk=>blk%blk_next
            end do
            blk_processor_tail=>blk
            if (associated(blk%blk_next)) then
                blk=>blk%blk_next
            end if
        end do
        !receive the block data from another processor
        blk=>blk_processor_head
        do i=1,np_nblk(rank+1)
            blk_id=blk%blk_id
            call block_id_to_processor_rank(np_nblk_old,blk_id,former_rank)
            if (rank/=former_rank) then
                if (allocated(blk%w)) then
                    call receive_block(blk,former_rank)
                else
                    call allocate_block_heavy_data(blk)
                    call receive_block(blk,former_rank)
                end if
            end if
            blk=>blk%blk_next
        end do
    else if (nd==2) then
    end if
    blk=>blk_head
    do i=1,nblk_total
        blk%on_processor=.false.
        blk=>blk%blk_next
    end do
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        blk%on_processor=.true.
        blk=>blk%blk_next
    end do
    deallocate(np_nblk_old)
end subroutine restart_smr_load_balancer

subroutine block_id_to_processor_rank(np_nblk_list,blk_id,rank_id)
    integer, dimension(:), allocatable :: np_nblk_list
    integer :: blk_id,rank_id,i,istart,iend
    istart=1;iend=np_nblk_list(1)
    do i=1,np
        if (blk_id>=istart.and.blk_id<=iend) then
            rank_id=i-1
            exit
        end if
        if (i/=np) then
            istart=istart+np_nblk_list(i)
            iend=iend+np_nblk_list(i+1)
        end if
    end do
end subroutine block_id_to_processor_rank

subroutine send_block(blk,receiver_rank)
    type(blockdef), pointer :: blk
    integer :: receiver_rank,tag,ierr,req
    tag=blk%blk_id*10+1
    call mpi_isend(blk%w,blk_primitive_cell_size,MPI_REAL8,receiver_rank,tag,MPI_COMM_WORLD,req,ierr)
    tag=blk%blk_id*10+2
    call mpi_isend(blk%u,blk_primitive_cell_size,MPI_REAL8,receiver_rank,tag,MPI_COMM_WORLD,req,ierr)
    tag=blk%blk_id*10+3
    call mpi_isend(blk%temp,blk_cell_size,MPI_REAL8,receiver_rank,tag,MPI_COMM_WORLD,req,ierr)
    tag=blk%blk_id*10+4
    call mpi_isend(blk%egv,blk_cell_size,MPI_REAL8,receiver_rank,tag,MPI_COMM_WORLD,req,ierr)
    if (iradiation/=0) then
        tag=blk%blk_id*10+5
        call mpi_isend(blk%Erad,blk_cell_size,MPI_REAL8,receiver_rank,tag,MPI_COMM_WORLD,req,ierr)
    end if
end subroutine send_block

subroutine receive_block(blk,former_rank)
    type(blockdef), pointer :: blk
    integer :: former_rank,tag,ierr,stat(MPI_STATUS_SIZE)
    tag=blk%blk_id*10+1
    call mpi_recv(blk%w,blk_primitive_cell_size,MPI_REAL8,former_rank,tag,MPI_COMM_WORLD,stat,ierr)
    tag=blk%blk_id*10+2
    call mpi_recv(blk%u,blk_primitive_cell_size,MPI_REAL8,former_rank,tag,MPI_COMM_WORLD,stat,ierr)
    tag=blk%blk_id*10+3
    call mpi_recv(blk%temp,blk_cell_size,MPI_REAL8,former_rank,tag,MPI_COMM_WORLD,stat,ierr)
    tag=blk%blk_id*10+4
    call mpi_recv(blk%egv,blk_cell_size,MPI_REAL8,former_rank,tag,MPI_COMM_WORLD,stat,ierr)
    if (iradiation/=0) then
        tag=blk%blk_id*10+5
        call mpi_recv(blk%Erad,blk_cell_size,MPI_REAL8,former_rank,tag,MPI_COMM_WORLD,stat,ierr)
    end if
end subroutine receive_block

subroutine grow_amr_node(blk)
    type(blockdef), pointer :: blk
    real(8) :: coords(2),n_size(2),zone(4)
    integer :: i,level
    if (nd==1) then
        if (refinement_criterion(blk)) then
            blk%derefine=.false.
            blk%discontinuity=.true.
            !recursively refine blk
            if (blk%level<max_refine_level) then
                call grow_tree_node_1d(blk)
                call record_processor_amr_grow_order(blk)
                blk%blk_xl%derefine=.false.
                blk%blk_xu%derefine=.false.
                call allocate_block_heavy_data(blk%blk_xl)
                call allocate_block_heavy_data(blk%blk_xu)
                call block_prolongation(blk)
                call grow_amr_node(blk%blk_xl)
                call grow_amr_node(blk%blk_xu)
            end if
            !refine the neighbouring blocks if a block satisfies the refinement criterion
            !if (associated(blk%blk_pre)) then
            !    blk_pre=>blk%blk_pre
            !    blk_pre%derefine=.false.
            !    if (blk_pre%level<max_refine_level) then
            !    end if
            !end if
            !if (associated(blk%blk_next)) then
            !    blk_next=>blk%blk_next
            !    blk%blk_next%derefine=.false.
            !end if
        else
            blk%derefine=derefinement_criterion(blk)
            blk%discontinuity=.false.
        end if
    else if (nd==2) then
    end if
end subroutine grow_amr_node

subroutine add_amr_buffer_blocks(blk)
    !see if blk's neighbours has more then 2 levels higher refinement
    !if so, refine this blk
    type(blockdef), pointer :: blk,blk_temp
    integer :: level
    if (nd==1) then
        level=blk%level
        if (associated(blk%blk_next)) then
            blk_temp=>blk%blk_next
            if (blk_temp%level>level+1) then
                call grow_amr_tree_node_buffer(blk)
                goto 10
            end if
        end if
        if (associated(blk%blk_pre)) then
            blk_temp=>blk%blk_pre
            if (blk_temp%level>level+1) then
                call grow_amr_tree_node_buffer(blk)
            end if
        end if
    else if (nd==2) then
10  end if
end subroutine add_amr_buffer_blocks

subroutine grow_amr_tree_node_buffer(blk)
    type(blockdef), pointer :: blk
    real(8) :: coords(2),n_size(2)
    integer :: id,level
    if (nd==1) then
        call grow_tree_node_1d(blk)
        call record_processor_amr_grow_order(blk)
        blk%derefine=.false.
        blk%blk_xl%derefine=.false.
        blk%blk_xu%derefine=.false.
        call allocate_block_heavy_data(blk%blk_xl)
        call allocate_block_heavy_data(blk%blk_xu)
        call block_prolongation(blk)
    else if (nd==2) then
    end if
end subroutine grow_amr_tree_node_buffer

subroutine derefine_amr_node(blk)
    !derefine the a block only when all its cousins are on this processor
    !and all of its cousins need to be derefined
    type(blockdef), pointer :: blk,blk_p
    logical :: l1,l2
    if (nd==1) then
        blk_p=>blk%blk_p
        if (associated(blk_p%blk_xl,blk_processor_tail)) then
            l1=.false.
        else
            l1=.true.
        end if
        if (associated(blk_p%blk_xu,blk_processor_head)) then
            l2=.false.
        else
            l2=.true.
        end if
        if (blk_p%blk_xl%derefine.and.blk_p%blk_xu%derefine.and.l1.and.l2) then   !derefine
            if (.not.allocated(blk_p%w)) then
                call allocate_block_heavy_data(blk_p)
            end if
            call block_restriction(blk_p)
            call record_processor_amr_derefine_order(blk_p)
            if (associated(blk_p%blk_xl,blk_processor_head)) then
                blk_processor_head=>blk_p
            end if
            if (associated(blk_p%blk_xu,blk_processor_tail)) then
                blk_processor_tail=>blk_p
            end if
            if (associated(blk%blk_pre)) then
                blk%blk_pre%blk_next=>blk_p
                blk_p%blk_pre=>blk%blk_pre
            else            !this is the first block
                blk_head=>blk_p
                blk_p%blk_next=>blk_p%blk_xu%blk_next
                blk_p%blk_xu%blk_next%blk_pre=>blk_p
            end if
            if (associated(blk%blk_next%blk_next)) then
                blk_p%blk_next=>blk_p%blk_xu%blk_next
                blk_p%blk_xu%blk_next%blk_pre=>blk_p
            else            !this is the last block
                blk_tail=>blk_p
                blk_p%blk_pre=>blk_p%blk_xl%blk_pre
                blk_p%blk_xl%blk_pre%blk_next=>blk_p
            end if
            call destroy_block(blk)
            call destroy_block(blk_p%blk_xu)
            nullify(blk_p%blk_xl,blk_p%blk_xu)
            blk=>blk_p
        end if
    else if (nd==2) then
    end if
end subroutine derefine_amr_node

function block_ownership(blk)
    !not working
    type(blockdef), pointer :: blk
    logical :: block_ownership
end function block_ownership

subroutine derefine_amr_node_conditional(blk)
    type(blockdef), pointer :: blk,blk_xl,blk_xu
    if (nd==1) then
        if (associated(blk%blk_xl)) then
            blk_xl=>blk%blk_xl;blk_xu=>blk%blk_xu
            if (associated(blk_xl,blk_processor_head)) then
                blk_processor_head=>blk
            end if
            if (associated(blk_xu,blk_processor_tail)) then
                blk_processor_tail=>blk
            end if
            if (associated(blk_xl%blk_pre)) then
                blk_xl%blk_pre%blk_next=>blk
                blk%blk_pre=>blk_xl%blk_pre
            else            !this is the first block
                blk_head=>blk
                blk%blk_next=>blk_xu%blk_next
                blk_xu%blk_next%blk_pre=>blk
            end if
            if (associated(blk_xu%blk_next)) then
                blk%blk_next=>blk_xu%blk_next
                blk_xu%blk_next%blk_pre=>blk
            else            !this is the last block
                blk_tail=>blk
                blk%blk_pre=>blk_xl%blk_pre
                blk_xl%blk_pre%blk_next=>blk
            end if
            call destroy_block(blk_xl)
            call destroy_block(blk_xu)
            nullify(blk%blk_xl,blk%blk_xu)
        end if
    else if (nd==2) then
    end if
end subroutine derefine_amr_node_conditional

function refinement_criterion(blk)
    type(blockdef), pointer :: blk
    logical :: refinement_criterion,shock,density,erad,temp
    real(8) :: pl,pr,vxl,vxr,egv,Etotal,gradp,trad1,trad2,trad3,xl,xc,xr
    integer :: i,j,k
    shock=.false.;density=.false.;erad=.false.;temp=.false.
    if (nd==1) then
        do i=0,blk_size_nx+1
            shock=shock_refine(blk%w(i-1,1,1,5),blk%w(i+1,1,1,5),blk%egv(i,1,1),blk%u(i,1,1,5))
            density=slope_refine(blk%w(i-1,1,1,1),blk%w(i,1,1,1),blk%w(i+1,1,1,1),rho_slope_refine_thresh)
            if (iradiation==4) then
                trad1=blk%temp(i-1,1,1)
                trad2=blk%temp(i,1,1)
                trad3=blk%temp(i+1,1,1)
                if (trad2>1.5d3) then
                    erad=slope_refine(trad1,trad2,trad3,erad_slope_refine_thresh)
                end if
            end if
            if (shock.or.density.or.erad.or.temp) then
                exit
            end if
        end do
        if (shock.or.density.or.erad.or.temp) then
            refinement_criterion=.true.
        else
            refinement_criterion=.false.
        end if
    else if (nd==2) then
    end if
end function refinement_criterion

function shock_refine(pl,pr,egv,Etotal)
    real(8) :: pl,pr,egv,Etotal,gradp
    logical :: shock_refine
    gradp=abs(pr-pl)/min(pl,pr)
    if (gradp>0.33.and.egv/Etotal>0.1) then
        shock_refine=.true.
    else
        shock_refine=.false.
    end if
end function shock_refine

function slope_refine(ql,qc,qr,slope_refine_thresh)
    real(8) :: ql,qc,qr,slope_refine_thresh
    logical :: slope_refine
    if (abs(qr-ql)/2/qc>slope_refine_thresh) then
        slope_refine=.true.
    else
        slope_refine=.false.
    end if
end function slope_refine

function second_order_gradient_refine()
    logical :: second_order_gradient_refine
end function second_order_gradient_refine

function derefinement_criterion(blk)
    !label each block's derefinement state
    type(blockdef), pointer :: blk
    real(8) :: pl,pr,vxl,vxr,egv,Etotal,gradp
    integer :: level,i,j,k
    integer, allocatable :: nb_level(:)
    logical :: derefinement_criterion,discon1,discon2
    if (blk%level>blk%static_level) then
        derefinement_criterion=.true.
        if (nd==1) then
            level=blk%level
            allocate(nb_level(2))
            nb_level=0
            if (associated(blk%blk_pre)) then
                if (associated(blk%blk_pre%blk_pre)) then
                    !nb_level(1)=max(blk%blk_pre%blk_pre%level,blk%blk_pre%level)
                    nb_level(1)=blk%blk_pre%level
                else
                    nb_level(1)=blk%blk_pre%level
                end if
                discon1=blk%blk_pre%discontinuity
            end if
            if (associated(blk%blk_next)) then
                if (associated(blk%blk_next%blk_next)) then
                    !nb_level(2)=max(blk%blk_next%blk_next%level,blk%blk_next%level)
                    nb_level(2)=blk%blk_next%level
                else
                    nb_level(2)=blk%blk_next%level
                end if
                discon2=blk%blk_next%discontinuity
            end if
            if (nb_level(1)-level>0.or.nb_level(2)-level>0.or.discon1.or.discon2) then
                !level difference cannot be greater than 1 between neighbouring blocks
                derefinement_criterion=.false.
            else
                !print *,level,nb_level,blk%key
                derefinement_criterion=.true.
            end if
            deallocate(nb_level)
        else if (nd==2) then
        end if
    else
        derefinement_criterion=.false.
    end if
end function derefinement_criterion

subroutine smr_load_balancer()
    !find the head and tail pointer on each processor
    !label the blocks that are on the processor
    type(blockdef), pointer :: blk
    integer :: blks_per_core,remainder,i,j
    call renumber_domain_blocks()
    blks_per_core=nblk_total/np
    np_nblk=blks_per_core
    remainder=mod(nblk_total,np)
    do i=1,remainder
        np_nblk(i)=np_nblk(i)+1
    end do
    blk=>blk_head
    do i=1,rank+1
        blk_processor_head=>blk
        do j=1,np_nblk(i)-1
            blk=>blk%blk_next
        end do
        blk_processor_tail=>blk
        if (associated(blk%blk_next)) then
            blk=>blk%blk_next
        end if
    end do
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        blk%on_processor=.true.
        blk=>blk%blk_next
    end do
end subroutine smr_load_balancer

subroutine examine_derefine()
    type(blockdef), pointer :: blk
    integer :: i
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        if (blk%derefine.eqv..true.) then
            blk%derefine=derefinement_criterion(blk)
        end if
        blk=>blk%blk_next
    end do
end subroutine examine_derefine

subroutine start_mpi()
    integer :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    if (ierr .ne. 0) then
        print*,'Unable to initialize PETSc'
        stop
    end if
end subroutine start_mpi

subroutine finalize_mpi()
    integer :: ierr
    call PetscFinalize(ierr)
end subroutine finalize_mpi

subroutine communicate_passive(sub)
    !passive scalars are associated with the matter
    procedure(extractarray) :: sub
    if (nd==1) then
    else if (nd==2) then
        call communicate_passive_internal(sub)
        call communicate_passive_external(sub)
    end if
end subroutine communicate_passive

subroutine communicate_passive_internal(sub)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    procedure(extractarray) :: sub
    real(8), dimension(:,:,:), pointer :: ap,ap_temp
    integer :: i,j,neighbour_rank
    blk=>blk_processor_head
    !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1),ap_temp(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    do i=1,np_nblk(rank+1)
        select case (blk%loc_type)
        case (1)
            call communicate_passive_edge_internal(blk,east,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,south,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,se,sub,ap,ap_temp)
        case (2)
            call communicate_passive_edge_internal(blk,east,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,west,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,south,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,se,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,sw,sub,ap,ap_temp)
        case (3)
            call communicate_passive_edge_internal(blk,west,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,south,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,sw,sub,ap,ap_temp)
        case (4)
            call communicate_passive_edge_internal(blk,east,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,south,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,north,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,se,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,ne,sub,ap,ap_temp)
        case (5)
            call communicate_passive_edge_internal(blk,east,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,west,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,south,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,north,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,se,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,sw,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,nw,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,ne,sub,ap,ap_temp)
        case (6)
            call communicate_passive_edge_internal(blk,west,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,south,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,north,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,sw,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,nw,sub,ap,ap_temp)
        case (7)
            call communicate_passive_edge_internal(blk,east,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,north,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,ne,sub,ap,ap_temp)
        case (8)
            call communicate_passive_edge_internal(blk,east,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,west,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,north,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,nw,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,ne,sub,ap,ap_temp)
        case (9)
            call communicate_passive_edge_internal(blk,west,sub,ap,ap_temp)
            call communicate_passive_edge_internal(blk,north,sub,ap,ap_temp)
            call communicate_passive_corner_internal(blk,nw,sub,ap,ap_temp)
        case default
        end select
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    nullify(ap,ap_temp)
end subroutine communicate_passive_internal

subroutine communicate_passive_edge_internal(blk,direction,sub,ap,ap_temp)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    procedure(extractarray) :: sub
    real(8), dimension(:,:,:), pointer :: ap,ap_temp
    integer :: neighbour_rank
    character(len=16) :: direction
    call link_a_neighbour_2d(blk,direction,blk_nb)
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank==rank) then
        blk_temp=>blk_nb%blk
        call sub(blk,ap)
        call sub(blk_temp,ap_temp)
        if (direction=='east') then
            ap(blk_size_nx+1:blk_size_nx+2,1:blk_size_ny,1)=ap_temp(1:2,1:blk_size_ny,1)
        else if (direction=='south') then
            ap(1:blk_size_nx,-1:0,1)=ap_temp(1:blk_size_nx,blk_size_ny-1:blk_size_ny,1)
        else if (direction=='west') then
            ap(-1:0,1:blk_size_ny,1)=ap_temp(blk_size_nx-1:blk_size_nx,1:blk_size_ny,1)
        else if (direction=='north') then
            ap(1:blk_size_nx,blk_size_ny+1:blk_size_ny+2,1)=ap_temp(1:blk_size_nx,1:2,1)
        end if
        nullify(blk_temp)
    end if
    nullify(blk_nb)
end subroutine communicate_passive_edge_internal

subroutine communicate_passive_corner_internal(blk,direction,sub,ap,ap_temp,mirror)
    type(blockdef), pointer :: blk,blk_temp
    type(blockneighbour), pointer :: blk_nb
    procedure(extractarray) :: sub
    real(8), dimension(:,:,:), pointer :: ap,ap_temp
    real(8), dimension(2,2) :: scalar_in,scalar_out
    integer :: neighbour_rank
    character(len=16) :: direction
    character(len=1), optional :: mirror
    logical :: transform
    if (present(mirror)) then
        transform=.true.
    else
        transform=.false.
    end if
    call link_a_neighbour_2d(blk,direction,blk_nb)
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank==rank) then
        blk_temp=>blk_nb%blk
        call sub(blk,ap)
        call sub(blk_temp,ap_temp)
        if (direction=='ne') then
            scalar_in=ap_temp(1:2,1:2,1)
        else if (direction=='se') then
            scalar_in=ap_temp(1:2,blk_size_ny-1:blk_size_ny,1)
        else if (direction=='sw') then
            scalar_in=ap_temp(blk_size_nx-1:blk_size_nx,blk_size_ny-1:blk_size_ny,1)
        else if (direction=='nw') then
            scalar_in=ap_temp(blk_size_nx-1:blk_size_nx,1:2,1)
        end if
        if (transform) then
            call corner_mirror_scalar(scalar_in,mirror,scalar_out)
        else
            scalar_out=scalar_in
        end if
        if (direction=='ne') then
            ap(blk_size_nx+1:blk_size_nx+2,blk_size_ny+1:blk_size_ny+2,1)=scalar_out
        else if (direction=='se') then
            ap(blk_size_nx+1:blk_size_nx+2,-1:0,1)=scalar_out
        else if (direction=='sw') then
            ap(-1:0,-1:0,1)=scalar_out
        else if (direction=='nw') then
            ap(-1:0,blk_size_ny+1:blk_size_ny+2,1)=scalar_out
        end if
        nullify(blk_temp)
    end if
    nullify(blk_nb)
end subroutine communicate_passive_corner_internal

subroutine communicate_passive_external(sub)
    type(blockdef), pointer :: blk
    procedure(extractarray) :: sub
    real(8), dimension(:,:,:), pointer :: ap
    integer :: i,ierr
    blk=>blk_processor_head
    !allocate(ap(blk_xlb:blk_xub,blk_ylb:blk_yub,1))
    do i=1,np_nblk(rank+1)
        call sub(blk,ap)
        select case (blk%loc_type)
        case (1)
            call communicate_passive_edge_external_send(blk,ap,east)
            call communicate_passive_edge_external_send(blk,ap,south)
            call communicate_passive_corner_external_send(blk,ap,se)
        case (2)
            call communicate_passive_edge_external_send(blk,ap,east)
            call communicate_passive_edge_external_send(blk,ap,west)
            call communicate_passive_edge_external_send(blk,ap,south)
            call communicate_passive_corner_external_send(blk,ap,se)
            call communicate_passive_corner_external_send(blk,ap,sw)
        case (3)
            call communicate_passive_edge_external_send(blk,ap,west)
            call communicate_passive_edge_external_send(blk,ap,south)
            call communicate_passive_corner_external_send(blk,ap,sw)
        case (4)
            call communicate_passive_edge_external_send(blk,ap,east)
            call communicate_passive_edge_external_send(blk,ap,south)
            call communicate_passive_edge_external_send(blk,ap,north)
            call communicate_passive_corner_external_send(blk,ap,se)
            call communicate_passive_corner_external_send(blk,ap,ne)
        case (5)
            call communicate_passive_edge_external_send(blk,ap,east)
            call communicate_passive_edge_external_send(blk,ap,west)
            call communicate_passive_edge_external_send(blk,ap,south)
            call communicate_passive_edge_external_send(blk,ap,north)
            call communicate_passive_corner_external_send(blk,ap,se)
            call communicate_passive_corner_external_send(blk,ap,sw)
            call communicate_passive_corner_external_send(blk,ap,nw)
            call communicate_passive_corner_external_send(blk,ap,ne)
        case (6)
            call communicate_passive_edge_external_send(blk,ap,west)
            call communicate_passive_edge_external_send(blk,ap,south)
            call communicate_passive_edge_external_send(blk,ap,north)
            call communicate_passive_corner_external_send(blk,ap,sw)
            call communicate_passive_corner_external_send(blk,ap,nw)
        case (7)
            call communicate_passive_edge_external_send(blk,ap,east)
            call communicate_passive_edge_external_send(blk,ap,north)
            call communicate_passive_corner_external_send(blk,ap,ne)
        case (8)
            call communicate_passive_edge_external_send(blk,ap,east)
            call communicate_passive_edge_external_send(blk,ap,west)
            call communicate_passive_edge_external_send(blk,ap,north)
            call communicate_passive_corner_external_send(blk,ap,nw)
            call communicate_passive_corner_external_send(blk,ap,ne)
        case (9)
            call communicate_passive_edge_external_send(blk,ap,west)
            call communicate_passive_edge_external_send(blk,ap,north)
            call communicate_passive_corner_external_send(blk,ap,nw)
        case default
        end select
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        call sub(blk,ap)
        select case (blk%loc_type)
        case (1)
            call communicate_passive_edge_external_recv(blk,ap,east)
            call communicate_passive_edge_external_recv(blk,ap,south)
            call communicate_passive_corner_external_recv(blk,ap,se)
        case (2)
            call communicate_passive_edge_external_recv(blk,ap,east)
            call communicate_passive_edge_external_recv(blk,ap,west)
            call communicate_passive_edge_external_recv(blk,ap,south)
            call communicate_passive_corner_external_recv(blk,ap,se)
            call communicate_passive_corner_external_recv(blk,ap,sw)
        case (3)
            call communicate_passive_edge_external_recv(blk,ap,west)
            call communicate_passive_edge_external_recv(blk,ap,south)
            call communicate_passive_corner_external_recv(blk,ap,sw)
        case (4)
            call communicate_passive_edge_external_recv(blk,ap,east)
            call communicate_passive_edge_external_recv(blk,ap,south)
            call communicate_passive_edge_external_recv(blk,ap,north)
            call communicate_passive_corner_external_recv(blk,ap,se)
            call communicate_passive_corner_external_recv(blk,ap,ne)
        case (5)
            call communicate_passive_edge_external_recv(blk,ap,east)
            call communicate_passive_edge_external_recv(blk,ap,west)
            call communicate_passive_edge_external_recv(blk,ap,south)
            call communicate_passive_edge_external_recv(blk,ap,north)
            call communicate_passive_corner_external_recv(blk,ap,se)
            call communicate_passive_corner_external_recv(blk,ap,sw)
            call communicate_passive_corner_external_recv(blk,ap,nw)
            call communicate_passive_corner_external_recv(blk,ap,ne)
        case (6)
            call communicate_passive_edge_external_recv(blk,ap,west)
            call communicate_passive_edge_external_recv(blk,ap,south)
            call communicate_passive_edge_external_recv(blk,ap,north)
            call communicate_passive_corner_external_recv(blk,ap,sw)
            call communicate_passive_corner_external_recv(blk,ap,nw)
        case (7)
            call communicate_passive_edge_external_recv(blk,ap,east)
            call communicate_passive_edge_external_recv(blk,ap,north)
            call communicate_passive_corner_external_recv(blk,ap,ne)
        case (8)
            call communicate_passive_edge_external_recv(blk,ap,east)
            call communicate_passive_edge_external_recv(blk,ap,west)
            call communicate_passive_edge_external_recv(blk,ap,north)
            call communicate_passive_corner_external_recv(blk,ap,nw)
            call communicate_passive_corner_external_recv(blk,ap,ne)
        case (9)
            call communicate_passive_edge_external_recv(blk,ap,west)
            call communicate_passive_edge_external_recv(blk,ap,north)
            call communicate_passive_corner_external_recv(blk,ap,nw)
        case default
        end select
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    nullify(ap)
end subroutine communicate_passive_external

subroutine communicate_passive_edge_external_send(blk,ap,direction)
    type(blockdef), pointer :: blk,blk_temp
    real(8), dimension(:,:,:), pointer :: ap
    real(8), dimension(:,:), allocatable :: msg_out
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    integer :: i,j,neighbour_rank,send_count,ierr,tag,stat(MPI_STATUS_SIZE),req
    call link_a_neighbour_2d(blk,direction,blk_nb)
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank/=rank) then
        blk_temp=>blk_nb%blk
        if (direction=='east') then
            send_count=2*blk_size_ny
            allocate(msg_out(2,blk_size_ny))
            msg_out(:,:)=ap(blk_size_nx-1:blk_size_nx,1:blk_size_ny,1)
        else if (direction=='south') then
            send_count=2*blk_size_nx
            allocate(msg_out(blk_size_nx,2))
            msg_out(:,:)=ap(1:blk_size_nx,1:2,1)
        else if (direction=='west') then
            send_count=2*blk_size_ny
            allocate(msg_out(2,blk_size_ny))
            msg_out(:,:)=ap(1:2,1:blk_size_ny,1)
        else if (direction=='north') then
            send_count=2*blk_size_nx
            allocate(msg_out(blk_size_nx,2))
            msg_out(:,:)=ap(1:blk_size_nx,blk_size_ny-1:blk_size_ny,1)
        end if
        tag=hash_tag_send(blk,blk_temp,direction)
        call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        call mpi_wait(req,stat,ierr)
        nullify(blk_temp)
        deallocate(msg_out)
    end if
end subroutine communicate_passive_edge_external_send

subroutine communicate_passive_corner_external_send(blk,ap,direction)
    type(blockdef), pointer :: blk,blk_temp
    real(8), dimension(:,:,:), pointer :: ap
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    integer :: i,j,neighbour_rank,send_count,ierr,tag,stat(MPI_STATUS_SIZE),req
    real(8) :: msg_out(2,2)
    call link_a_neighbour_2d(blk,direction,blk_nb)
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank/=rank) then
        blk_temp=>blk_nb%blk
        send_count=4
        if (direction=='se') then
            msg_out=ap(blk_size_nx-1:blk_size_nx,1:2,1)
        else if (direction=='sw') then
            msg_out=ap(1:2,1:2,1)
        else if (direction=='nw') then
            msg_out=ap(1:2,blk_size_ny-1:blk_size_ny,1)
        else if (direction=='ne') then
            msg_out=ap(blk_size_nx-1:blk_size_nx,blk_size_ny-1:blk_size_ny,1)
        end if
        tag=hash_tag_send(blk,blk_temp,direction)
        call mpi_isend(msg_out,send_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,req,ierr)
        call mpi_wait(req,stat,ierr)
        nullify(blk_temp)
    end if
end subroutine communicate_passive_corner_external_send

subroutine communicate_passive_edge_external_recv(blk,ap,direction)
    type(blockdef), pointer :: blk,blk_temp
    real(8), dimension(:,:,:), pointer :: ap
    real(8), dimension(:,:), allocatable :: msg_in
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    integer :: i,j,neighbour_rank,recv_count,ierr,tag,stat(MPI_STATUS_SIZE),req
    call link_a_neighbour_2d(blk,direction,blk_nb)
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank/=rank) then
        blk_temp=>blk_nb%blk
        if (direction=='east'.or.direction=='west') then
            recv_count=2*blk_size_ny
            allocate(msg_in(2,blk_size_ny))
        else if (direction=='south'.or.direction=='north') then
            recv_count=2*blk_size_nx
            allocate(msg_in(blk_size_nx,2))
        end if
        tag=hash_tag_recv(blk_temp,blk,direction)
        call mpi_recv(msg_in,recv_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        if (direction=='east') then
            ap(blk_size_nx+1:blk_size_nx+2,1:blk_size_ny,1)=msg_in(:,:)
        else if (direction=='south') then
            ap(1:blk_size_nx,-1:0,1)=msg_in(:,:)
        else if (direction=='west') then
            ap(-1:0,1:blk_size_ny,1)=msg_in(:,:)
        else if (direction=='north') then
            ap(1:blk_size_nx,blk_size_ny+1:blk_size_ny+2,1)=msg_in(:,:)
        end if
        nullify(blk_temp)
        deallocate(msg_in)
    end if
end subroutine communicate_passive_edge_external_recv

subroutine communicate_passive_corner_external_recv(blk,ap,direction,mirror)
    !scalar can have mirror operation, but not reflect
    type(blockdef), pointer :: blk,blk_temp
    real(8), dimension(:,:,:), pointer :: ap
    real(8) :: msg_in(2,2)
    type(blockneighbour), pointer :: blk_nb
    character(len=16) :: direction
    integer :: i,j,neighbour_rank,recv_count,ierr,tag,stat(MPI_STATUS_SIZE),req
    character(len=1), optional :: mirror
    call link_a_neighbour_2d(blk,direction,blk_nb)
    neighbour_rank=blk_nb%cpu_rank
    if (neighbour_rank/=rank) then
        blk_temp=>blk_nb%blk
        recv_count=4
        tag=hash_tag_recv(blk_temp,blk,direction)
        call mpi_recv(msg_in,recv_count,MPI_REAL8,neighbour_rank,tag,MPI_COMM_WORLD,stat,ierr)
        if (direction=='se') then
            ap(blk_size_nx+1:blk_size_nx+2,-1:0,1)=msg_in
        else if (direction=='sw') then
            ap(-1:0,-1:0,1)=msg_in
        else if (direction=='nw') then
            ap(-1:0,blk_size_ny+1:blk_size_ny+2,1)=msg_in
        else if (direction=='ne') then
            ap(blk_size_nx+1:blk_size_nx+2,blk_size_ny+1:blk_size_ny+2,1)=msg_in
        end if
        nullify(blk_temp)
    end if
end subroutine communicate_passive_corner_external_recv

end module communication
