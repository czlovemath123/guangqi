module petsc_fld_2d
#include <petsc/finclude/petscksp.h>
use datastructure
use eos
use communication
use phylib
use radiation_common_functions
use petscksp
implicit none

contains

subroutine petsc_assemble_mat_vec_2d(dt,A,b)
    type(blockdef), pointer :: blk_traversal
    Mat A
    Vec b
    real(8) :: dt
    integer :: i
    blk_traversal=>blk_processor_head
    do i=1,np_nblk(rank+1)
        select case (blk_traversal%loc_type)
        case (1)
            call petsc_Ab_loc1_2d_block(blk_traversal,dt,A,b)
        case (2)
            call petsc_Ab_loc2_2d_block(blk_traversal,dt,A,b)
        case (3)
            call petsc_Ab_loc3_2d_block(blk_traversal,dt,A,b)
        case (4)
            call petsc_Ab_loc4_2d_block(blk_traversal,dt,A,b)
        case (5)
            call petsc_Ab_loc5_2d_block(blk_traversal,dt,A,b)
        case (6)
            call petsc_Ab_loc6_2d_block(blk_traversal,dt,A,b)
        case (7)
            call petsc_Ab_loc7_2d_block(blk_traversal,dt,A,b)
        case (8)
            call petsc_Ab_loc8_2d_block(blk_traversal,dt,A,b)
        case (9)
            call petsc_Ab_loc9_2d_block(blk_traversal,dt,A,b)
        end select
        if (i/=np_nblk(rank+1)) then
            blk_traversal=>blk_traversal%blk_next
        end if
    end do
end subroutine petsc_assemble_mat_vec_2d

subroutine petsc_Ab_loc1_2d_block(blk,dt,A,b)
    !blk%loc_type=1
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3)
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    i=1;j=blk_size_ny
    call petsc_Ab_1_2d_cell(blk,dt,A,b,i,j)
    j=blk_size_ny
    do i=2,blk_size_nx
        call petsc_Ab_2_2d_cell(blk,dt,A,b,i,j)
    end do
    i=1
    do j=1,blk_size_ny-1
        call petsc_Ab_4_2d_cell(blk,dt,A,b,i,j)
    end do
    do j=1,blk_size_ny-1
        do i=2,blk_size_nx
            call petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
        end do
    end do
end subroutine petsc_Ab_loc1_2d_block

subroutine petsc_Ab_loc2_2d_block(blk,dt,A,b)
    !blk%loc_type=2
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3)
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    j=blk_size_ny
    do i=1,blk_size_nx
        call petsc_Ab_2_2d_cell(blk,dt,A,b,i,j)
    end do
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx
            call petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
        end do
    end do
end subroutine petsc_Ab_loc2_2d_block

subroutine petsc_Ab_loc3_2d_block(blk,dt,A,b)
    !blk%loc_type=3
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3)
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    i=blk_size_nx;j=blk_size_ny
    call petsc_Ab_3_2d_cell(blk,dt,A,b,i,j)
    j=blk_size_ny
    do i=1,blk_size_nx-1
        call petsc_Ab_2_2d_cell(blk,dt,A,b,i,j)
    end do
    i=blk_size_nx
    do j=1,blk_size_ny-1
        call petsc_Ab_6_2d_cell(blk,dt,A,b,i,j)
    end do
    do j=1,blk_size_ny-1
        do i=1,blk_size_nx-1
            call petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
        end do
    end do
end subroutine petsc_Ab_loc3_2d_block

subroutine petsc_Ab_loc4_2d_block(blk,dt,A,b)
    !blk%loc_type=4
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3)
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    i=1
    do j=1,blk_size_ny
        call petsc_Ab_4_2d_cell(blk,dt,A,b,i,j)
    end do
    do j=1,blk_size_ny
        do i=2,blk_size_nx
            call petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
        end do
    end do
end subroutine petsc_Ab_loc4_2d_block

subroutine petsc_Ab_loc5_2d_block(blk,dt,A,b)
    !blk%loc_type=5
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3)
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            call petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
        end do
    end do
end subroutine petsc_Ab_loc5_2d_block

subroutine petsc_Ab_loc6_2d_block(blk,dt,A,b)
    !blk%loc_type=6
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3)
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    i=blk_size_nx
    do j=1,blk_size_ny
        call petsc_Ab_6_2d_cell(blk,dt,A,b,i,j)
    end do
    do j=1,blk_size_ny
        do i=1,blk_size_nx-1
            call petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
        end do
    end do
end subroutine petsc_Ab_loc6_2d_block

subroutine petsc_Ab_loc7_2d_block(blk,dt,A,b)
    !blk%loc_type=7
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3)
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    i=1;j=1
    call petsc_Ab_7_2d_cell(blk,dt,A,b,i,j)
    i=1
    do j=2,blk_size_ny
        call petsc_Ab_4_2d_cell(blk,dt,A,b,i,j)
    end do
    j=1
    do i=2,blk_size_nx
        call petsc_Ab_8_2d_cell(blk,dt,A,b,i,j)
    end do
    do j=2,blk_size_ny
        do i=2,blk_size_nx
            call petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
        end do
    end do
end subroutine petsc_Ab_loc7_2d_block

subroutine petsc_Ab_loc8_2d_block(blk,dt,A,b)
    !blk%loc_type=8
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3)
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    j=1
    do i=1,blk_size_nx
        call petsc_Ab_8_2d_cell(blk,dt,A,b,i,j)
    end do
    do j=2,blk_size_ny
        do i=1,blk_size_nx
            call petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
        end do
    end do
end subroutine petsc_Ab_loc8_2d_block

subroutine petsc_Ab_loc9_2d_block(blk,dt,A,b)
    !blk%loc_type=9
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3)
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    i=blk_size_nx;j=1
    call petsc_Ab_9_2d_cell(blk,dt,A,b,i,j)
    i=blk_size_nx
    do j=2,blk_size_ny
        call petsc_Ab_6_2d_cell(blk,dt,A,b,i,j)
    end do
    j=1
    do i=1,blk_size_nx-1
        call petsc_Ab_8_2d_cell(blk,dt,A,b,i,j)
    end do
    do j=2,blk_size_ny
        do i=1,blk_size_nx-1
            call petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
        end do
    end do
end subroutine petsc_Ab_loc9_2d_block

subroutine petsc_Ab_1_2d_cell(blk,dt,A,b,i,j)
    !cell type is upper left corner
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3),Dx_l,Dx_u,Dy_l,Dy_u,temp,q
    integer :: i,j,ii_E,ii_T,ii,jj,ierr
    ii=blkid_ij_to_matrix_ii(blk,i,j)
    ii_E=ii
    ii_T=ii_E+1
    dx=blk%dxyz(1)
    alpha=dt/dx/dx
    beta=blk%sigma_planck(i,j,1)*c_light*dt
    temp=blk%temp(i,j,1)
    Dx_u=blk%kx(i,j,1)
    Dy_l=blk%ky(i,j-1,1)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4d0*beta*a_rad*temp**3+blk%cv(i,j,1)
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4d0*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i+1,j)
    v=-alpha*Dx_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j-1)
    v=-alpha*Dy_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (rad_bound_type(1)==2.and.rad_bound_type(4)==2) then
        v=1d0+(beta+alpha*(Dx_u+Dy_l))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else if (rad_bound_type(1)==1.and.rad_bound_type(4)==1) then
        if (igeometry==0) then
            Dx_l=blk%kx(i-1,j,1)
            Dy_u=blk%ky(i,j,1)
            q=1.5d0*blk%sigma_rosseland(i,j,1)*dx
            v=1d0+(beta+alpha*(q*Dx_l/(1+q)+Dx_u+Dy_l+q*Dy_u/(1+q)))*const_Erad
            call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        else if (igeometry==1) then
        else if (igeometry==2) then
        end if
    else
        print *,'2d fld under construction'
        stop
    end if
    v=blk%Erad(i,j,1)-(3*beta*a_rad*blk%temp(i,j,1)**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%cv(i,j,1)*blk%temp(i,j,1)+3*beta*a_rad*blk%temp(i,j,1)**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_1_2d_cell

subroutine petsc_Ab_2_2d_cell(blk,dt,A,b,i,j)
    !cell type is upper edge
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3),Dx_l,Dx_u,Dy_l,Dy_u,temp,q
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    ii=blkid_ij_to_matrix_ii(blk,i,j)
    ii_E=ii
    ii_T=ii_E+1
    dx=blk%dxyz(1)
    alpha=dt/dx/dx
    beta=blk%sigma_planck(i,j,1)*c_light*dt
    temp=blk%temp(i,j,1)
    Dx_l=blk%kx(i-1,j,1)
    Dx_u=blk%kx(i,j,1)
    Dy_l=blk%ky(i,j-1,1)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4d0*beta*a_rad*temp**3+blk%cv(i,j,1)
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4d0*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i-1,j)
    v=-alpha*Dx_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i+1,j)
    v=-alpha*Dx_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j-1)
    v=-alpha*Dy_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (rad_bound_type(4)==2) then
        v=1d0+(beta+alpha*(Dx_l+Dx_u+Dy_l))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else if (rad_bound_type(4)==1) then
        if (igeometry==0) then
            Dy_u=blk%ky(i,j,1)
            q=1.5d0*blk%sigma_rosseland(i,j,1)*dx
            v=1d0+(beta+alpha*(Dx_l+Dx_u+Dy_l+q*Dy_u/(1+q)))*const_Erad
            call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        else if (igeometry==1) then
        else if (igeometry==2) then
        end if
    else
        print *,'2d fld under construction'
        stop
    end if
    v=blk%Erad(i,j,1)-(3*beta*a_rad*blk%temp(i,j,1)**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%cv(i,j,1)*blk%temp(i,j,1)+3*beta*a_rad*blk%temp(i,j,1)**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_2_2d_cell

subroutine petsc_Ab_3_2d_cell(blk,dt,A,b,i,j)
    !cell type is upper right corner
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3),Dx_l,Dx_u,Dy_l,Dy_u,temp,q
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    ii=blkid_ij_to_matrix_ii(blk,i,j)
    ii_E=ii
    ii_T=ii_E+1
    dx=blk%dxyz(1)
    alpha=dt/dx/dx
    beta=blk%sigma_planck(i,j,1)*c_light*dt
    temp=blk%temp(i,j,1)
    Dx_l=blk%kx(i-1,j,1)
    Dy_l=blk%ky(i,j-1,1)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4d0*beta*a_rad*temp**3+blk%cv(i,j,1)
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4d0*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i-1,j)
    v=-alpha*Dx_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j-1)
    v=-alpha*Dy_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (rad_bound_type(2)==2.and.rad_bound_type(4)==2) then
        v=1d0+(beta+alpha*(Dx_l+Dy_l))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else if (rad_bound_type(2)==1.and.rad_bound_type(4)==1) then
        if (igeometry==0) then
            Dx_u=blk%kx(i,j,1)
            Dy_u=blk%ky(i,j,1)
            q=1.5d0*blk%sigma_rosseland(i,j,1)*dx
            v=1d0+(beta+alpha*(Dx_l+q*Dx_u/(1+q)+Dy_l+q*Dy_u/(1+q)))*const_Erad
            call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        else if (igeometry==1) then
        else if (igeometry==2) then
        end if
    else
        print *,'2d fld under construction'
        stop
    end if
    v=blk%Erad(i,j,1)-(3*beta*a_rad*blk%temp(i,j,1)**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%cv(i,j,1)*blk%temp(i,j,1)+3*beta*a_rad*blk%temp(i,j,1)**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_3_2d_cell

subroutine petsc_Ab_4_2d_cell(blk,dt,A,b,i,j)
    !cell type is left edge
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3),Dx_l,Dx_u,Dy_l,Dy_u,temp,q
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    ii=blkid_ij_to_matrix_ii(blk,i,j)
    ii_E=ii
    ii_T=ii_E+1
    dx=blk%dxyz(1)
    alpha=dt/dx/dx
    beta=blk%sigma_planck(i,j,1)*c_light*dt
    temp=blk%temp(i,j,1)
    Dx_u=blk%kx(i,j,1)
    Dy_l=blk%ky(i,j-1,1)
    Dy_u=blk%ky(i,j,1)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4d0*beta*a_rad*temp**3+blk%cv(i,j,1)
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4d0*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i+1,j)
    v=-alpha*Dx_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j-1)
    v=-alpha*Dy_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j+1)
    v=-alpha*Dy_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (rad_bound_type(1)==2) then
        v=1d0+(beta+alpha*(Dx_u+Dy_l+Dy_u))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else if (rad_bound_type(1)==1) then
        if (igeometry==0) then
            Dx_l=blk%kx(i-1,j,1)
            q=1.5d0*blk%sigma_rosseland(i,j,1)*dx
            v=1d0+(beta+alpha*(q*Dx_l/(1+q)+Dx_u+Dy_l+Dy_u))*const_Erad
            call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        else if (igeometry==1) then
        else if (igeometry==2) then
        end if
    else
        print *,'2d fld under construction'
        stop
    end if
    v=blk%Erad(i,j,1)-(3*beta*a_rad*blk%temp(i,j,1)**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%cv(i,j,1)*blk%temp(i,j,1)+3*beta*a_rad*blk%temp(i,j,1)**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_4_2d_cell

subroutine petsc_Ab_5_2d_cell(blk,dt,A,b,i,j)
    !cell type is interior
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3),Dx_l,Dx_u,Dy_l,Dy_u,temp,q
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    ii=blkid_ij_to_matrix_ii(blk,i,j)
    ii_E=ii
    ii_T=ii_E+1
    dx=blk%dxyz(1)
    alpha=dt/dx/dx
    beta=blk%sigma_planck(i,j,1)*c_light*dt
    temp=blk%temp(i,j,1)
    Dx_l=blk%kx(i-1,j,1)
    Dx_u=blk%kx(i,j,1)
    Dy_l=blk%ky(i,j-1,1)
    Dy_u=blk%ky(i,j,1)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4d0*beta*a_rad*temp**3+blk%cv(i,j,1)
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4d0*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i-1,j)
    v=-alpha*Dx_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i+1,j)
    v=-alpha*Dx_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j-1)
    v=-alpha*Dy_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j+1)
    v=-alpha*Dy_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=1d0+(beta+alpha*(Dx_l+Dx_u+Dy_l+Dy_u))*const_Erad
    call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%Erad(i,j,1)-(3*beta*a_rad*blk%temp(i,j,1)**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%cv(i,j,1)*blk%temp(i,j,1)+3*beta*a_rad*blk%temp(i,j,1)**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_5_2d_cell

subroutine petsc_Ab_6_2d_cell(blk,dt,A,b,i,j)
    !cell type is right edge
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3),Dx_l,Dx_u,Dy_l,Dy_u,temp,q
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    ii=blkid_ij_to_matrix_ii(blk,i,j)
    ii_E=ii
    ii_T=ii_E+1
    dx=blk%dxyz(1)
    alpha=dt/dx/dx
    beta=blk%sigma_planck(i,j,1)*c_light*dt
    temp=blk%temp(i,j,1)
    Dx_l=blk%kx(i-1,j,1)
    Dy_l=blk%ky(i,j-1,1)
    Dy_u=blk%ky(i,j,1)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4d0*beta*a_rad*temp**3+blk%cv(i,j,1)
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4d0*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i-1,j)
    v=-alpha*Dx_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j-1)
    v=-alpha*Dy_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j+1)
    v=-alpha*Dy_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (rad_bound_type(2)==2) then
        v=1d0+(beta+alpha*(Dx_l+Dy_l+Dy_u))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else if (rad_bound_type(2)==1) then
        if (igeometry==0) then
            Dx_u=blk%kx(i,j,1)
            q=1.5d0*blk%sigma_rosseland(i,j,1)*dx
            v=1d0+(beta+alpha*(Dx_l+q*Dx_u/(1+q)+Dy_l+Dy_u))*const_Erad
            call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
        else if (igeometry==1) then
        else if (igeometry==2) then
        end if
    else
        print *,'2d fld under construction'
        stop
    end if
    v=blk%Erad(i,j,1)-(3*beta*a_rad*blk%temp(i,j,1)**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%cv(i,j,1)*blk%temp(i,j,1)+3*beta*a_rad*blk%temp(i,j,1)**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_6_2d_cell

subroutine petsc_Ab_7_2d_cell(blk,dt,A,b,i,j)
    !cell type is lower left corner
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3),Dx_l,Dx_u,Dy_l,Dy_u,temp,q
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    ii=blkid_ij_to_matrix_ii(blk,i,j)
    ii_E=ii
    ii_T=ii_E+1
    dx=blk%dxyz(1)
    alpha=dt/dx/dx
    beta=blk%sigma_planck(i,j,1)*c_light*dt
    temp=blk%temp(i,j,1)
    Dx_u=blk%kx(i,j,1)
    Dy_u=blk%ky(i,j,1)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4d0*beta*a_rad*temp**3+blk%cv(i,j,1)
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4d0*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i+1,j)
    v=-alpha*Dx_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j+1)
    v=-alpha*Dy_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (rad_bound_type(1)==2.and.rad_bound_type(3)==2) then
        v=1d0+(beta+alpha*(Dx_u+Dy_u))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else if (rad_bound_type(1)==1.and.rad_bound_type(3)==1) then
        Dx_l=blk%kx(i-1,j,1)
        Dy_l=blk%ky(i,j-1,1)
        q=1.5d0*blk%sigma_rosseland(i,j,1)*dx
        v=1d0+(beta+alpha*(q*Dx_l/(1+q)+Dx_u+q*Dy_l/(1+q)+Dy_u))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else
        print *,'2d fld under construction'
        stop
    end if
    v=blk%Erad(i,j,1)-(3*beta*a_rad*blk%temp(i,j,1)**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%cv(i,j,1)*blk%temp(i,j,1)+3*beta*a_rad*blk%temp(i,j,1)**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_7_2d_cell

subroutine petsc_Ab_8_2d_cell(blk,dt,A,b,i,j)
    !cell type is lower edge
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3),Dx_l,Dx_u,Dy_l,Dy_u,temp,q
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    ii=blkid_ij_to_matrix_ii(blk,i,j)
    ii_E=ii
    ii_T=ii_E+1
    dx=blk%dxyz(1)
    alpha=dt/dx/dx
    beta=blk%sigma_planck(i,j,1)*c_light*dt
    temp=blk%temp(i,j,1)
    Dx_l=blk%kx(i-1,j,1)
    Dx_u=blk%kx(i,j,1)
    Dy_u=blk%ky(i,j,1)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4d0*beta*a_rad*temp**3+blk%cv(i,j,1)
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4d0*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i-1,j)
    v=-alpha*Dx_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i+1,j)
    v=-alpha*Dx_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j+1)
    v=-alpha*Dy_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (rad_bound_type(3)==2) then
        v=1d0+(beta+alpha*(Dx_l+Dx_u+Dy_u))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else if (rad_bound_type(3)==1) then
        Dy_l=blk%ky(i,j-1,1)
        q=1.5d0*blk%sigma_rosseland(i,j,1)*dx
        v=1d0+(beta+alpha*(Dx_l+Dx_u+q*Dy_l/(1+q)+Dy_u))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else
        print *,'2d fld under construction'
        stop
    end if
    v=blk%Erad(i,j,1)-(3*beta*a_rad*blk%temp(i,j,1)**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%cv(i,j,1)*blk%temp(i,j,1)+3*beta*a_rad*blk%temp(i,j,1)**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_8_2d_cell

subroutine petsc_Ab_9_2d_cell(blk,dt,A,b,i,j)
    !cell type is lower right corner
    type(blockdef), pointer :: blk
    Mat A
    Vec b
    real(8) :: v,alpha,beta,dt,dx,blk_dxyz(3),Dx_l,Dx_u,Dy_l,Dy_u,temp,q
    integer :: i,j,ii_E,ii_T,ii,jj,ierr,blk_id
    ii=blkid_ij_to_matrix_ii(blk,i,j)
    ii_E=ii
    ii_T=ii_E+1
    dx=blk%dxyz(1)
    alpha=dt/dx/dx
    beta=blk%sigma_planck(i,j,1)*c_light*dt
    temp=blk%temp(i,j,1)
    Dx_l=blk%kx(i-1,j,1)
    Dy_u=blk%ky(i,j,1)
    v=-beta
    call MatSetValue(A,ii_T,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=4d0*beta*a_rad*temp**3+blk%cv(i,j,1)
    call MatSetValue(A,ii_T,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=-4d0*beta*a_rad*temp**3*const_Erad
    call MatSetValue(A,ii_E,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i-1,j)
    v=-alpha*Dx_l*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    jj=blkid_ij_to_matrix_ii(blk,i,j+1)
    v=-alpha*Dy_u*const_Erad
    call MatSetValue(A,ii_E,jj,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (rad_bound_type(2)==2.and.rad_bound_type(3)==2) then
        v=1d0+(beta+alpha*(Dx_l+Dy_u))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else if (rad_bound_type(2)==1.and.rad_bound_type(3)==1) then
        Dx_u=blk%kx(i,j,1)
        Dy_l=blk%ky(i,j-1,1)
        q=1.5d0*blk%sigma_rosseland(i,j,1)*dx
        v=1d0+(beta+alpha*(Dx_l+q*Dx_u/(1+q)+q*Dy_l/(1+q)+Dy_u))*const_Erad
        call MatSetValue(A,ii_E,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    else
        print *,'2d fld under construction'
        stop
    end if
    v=blk%Erad(i,j,1)-(3*beta*a_rad*blk%temp(i,j,1)**4)*const_Erad
    call VecSetValue(b,ii_E,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
    v=blk%cv(i,j,1)*blk%temp(i,j,1)+3*beta*a_rad*blk%temp(i,j,1)**4
    call VecSetValue(b,ii_T,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
end subroutine petsc_Ab_9_2d_cell

function check_interior_cell(blk,i,j)
    !true if the cell is not a domain boundary guard cell
    !true if the cell is a block guard cell but not a domain boundary guard cell
    type(blockdef), pointer :: blk
    integer :: i,j
    logical :: check_interior_cell
    if (blk%loc_type==5) then
        check_interior_cell=.true.
    else
        select case (blk%loc_type)
        case (1)
            if (i<=0.or.j>=blk_size_ny+1) then
                check_interior_cell=.false.
            else
                check_interior_cell=.true.
            end if
        case (2)
            if (j>=blk_size_ny+1) then
                check_interior_cell=.false.
            else
                check_interior_cell=.true.
            end if
        case (3)
            if (i>=blk_size_nx+1.or.j>=blk_size_ny+1) then
                check_interior_cell=.false.
            else
                check_interior_cell=.true.
            end if
        case (4)
            if (i<=0) then
                check_interior_cell=.false.
            else
                check_interior_cell=.true.
            end if
        case (6)
            if (i>=blk_size_nx+1) then
                check_interior_cell=.false.
            else
                check_interior_cell=.true.
            end if
        case (7)
            if (i<=0.or.j<=0) then
                check_interior_cell=.false.
            else
                check_interior_cell=.true.
            end if
        case (8)
            if (j<=0) then
                check_interior_cell=.false.
            else
                check_interior_cell=.true.
            end if
        case (9)
            if (i>=blk_size_nx+1.or.j<=0) then
                check_interior_cell=.false.
            else
                check_interior_cell=.true.
            end if
        case default
            print *,'not such loc_type'
            stop
        end select
    end if
end function check_interior_cell

function blkid_ij_to_matrix_ii(blk,i,j)
    !for interior cells only (not boundary cells)
    type(blockdef), pointer :: blk,blk_nb
    type(blockneighbour), pointer :: nb
    integer :: blk_id,i,j,blkid_ij_to_matrix_ii,i_nb,j_nb,blk_id_nb
    logical :: interior
    interior=check_interior_cell(blk,i,j)
    if (.not.interior) then
        print *,'not interior cells, check matrix assembling process'
        print *,blk%blk_id,blk%loc_type,i,j
        stop
    end if
    if (i>=1.and.i<=blk_size_nx.and.j>=1.and.j<=blk_size_ny) then
        blk_id=blk%blk_id
        blkid_ij_to_matrix_ii=((blk_id-1)*blk_size_nx*blk_size_ny+(j-1)*blk_size_nx+i-1)*2
    else
        if (i==0.and.j==0) then
            !sw corner
            nb=>blk%nb_sw
            blk_id_nb=nb%blk%blk_id
            i_nb=blk_size_nx
            j_nb=blk_size_ny
        else if (i==0.and.j>=1.and.j<=blk_size_ny) then
            !west edge
            nb=>blk%nb_w
            blk_id_nb=nb%blk%blk_id
            i_nb=blk_size_nx
            j_nb=j
        else if (i==0.and.j==blk_size_ny+1) then
            !nw corner
            nb=>blk%nb_nw
            blk_id_nb=nb%blk%blk_id
            i_nb=blk_size_nx
            j_nb=1
        else if (i==blk_size_nx+1.and.j==0) then
            !se corner
            nb=>blk%nb_se
            blk_id_nb=nb%blk%blk_id
            i_nb=1
            j_nb=blk_size_ny
        else if (i==blk_size_nx+1.and.j==blk_size_ny+1) then
            !ne corner
            nb=>blk%nb_ne
            blk_id_nb=nb%blk%blk_id
            i_nb=1
            j_nb=1
        else if (i==blk_size_nx+1.and.j>=1.and.j<=blk_size_ny) then
            !east edge
            nb=>blk%nb_e
            blk_id_nb=nb%blk%blk_id
            i_nb=1
            j_nb=j
        else if (j==0.and.i>=1.and.i<=blk_size_nx) then
            !south edge
            nb=>blk%nb_s
            blk_id_nb=nb%blk%blk_id
            i_nb=i
            j_nb=blk_size_ny
        else if (j==blk_size_ny+1.and.i>=1.and.i<=blk_size_nx) then
            !north edge
            nb=>blk%nb_n
            blk_id_nb=nb%blk%blk_id
            i_nb=i
            j_nb=1
        end if
        blkid_ij_to_matrix_ii=((blk_id_nb-1)*blk_size_nx*blk_size_ny+(j_nb-1)*blk_size_nx+i_nb-1)*2
    end if
end function blkid_ij_to_matrix_ii

subroutine matrix_ii_to_blkid_ij(ii,blk_id,i,j)
    integer :: ii,blk_id,i,j,blk_size,row_size,v1
    blk_size=2*blk_size_nx*blk_size_ny
    row_size=2*blk_size_nx
    if (mod(ii,blk_size)==0) then
        blk_id=ii/blk_size
        i=blk_size_nx
        j=blk_size_ny
    else
        blk_id=ii/blk_size+1
        v1=mod(ii,blk_size)
        if (mod(v1,row_size)==0) then
            j=v1/row_size
            i=blk_size_nx
        else
            j=v1/row_size+1
            i=mod(v1,row_size)/2
        end if
    end if
end subroutine matrix_ii_to_blkid_ij

end module petsc_fld_2d
