module godunov
use phylib
use hydro
use eos
use datastructure
use source_control
use radiation
implicit none

contains

function t_hydro_block(vglobalmax,dx)
    real(8) :: vglobalmax(nd),t_hydro_block,dx
    t_hydro_block=dx*CFL/maxval(vglobalmax)
    if (time_sys%t+t_hydro_block>=time_sys%t_next) then
        t_hydro_block=time_sys%t_next-time_sys%t
    end if
end function t_hydro_block

subroutine godunov_hydro_split()
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:,:), allocatable :: du,xfluxarray,yfluxarray,dsource
    real(8), dimension(:,:,:,:), allocatable :: du_x,du_y,dudt_x,dudt_y,u1,u2,w1,w2
    real(8), dimension(:,:,:), allocatable :: temp1,temp2,egv1,egv2,dqdt_x,dqdt_y,q1,q2
    real(8), dimension(5) :: wl,wr,xflux,yflux,s_geo1,s_geo2,source,u_temp
    real(8) :: vglobalmax(nd),dx
    real(8), dimension(:), allocatable :: dt_hydro_block
    integer :: i,j,k,l,m,n,left,right,ijk(3)
    character(len=1) :: dir
    character(len=128) :: alert
    real(8) :: temp(2),egv(2)
    allocate(dt_hydro_block(nblk_total))
    if (nd==1) then
        call allocate_xsurface_data_block(xfluxarray,5)
        call allocate_cell_data_block(du,5)
        blk=>blk_head
        do while (associated(blk))
            vglobalmax=0d0
            dx=blk%dxyz(1)
            if (igeometry/=0) then
                call geometry_source(blk,blk%s_geo1)
            end if
#if         isolver==1
            call hllc_sweep(blk%w,vglobalmax,dx,xfluxarray)
#elif       isolver==2
            call eos_hllc_analytic_sweep(blk%w,blk%temp,blk%egv,  &
                blk%dudt_x,vglobalmax,dx,xfluxarray)
#endif
            if (iradiation/=0) then
                blk%mass_flux_x=xfluxarray(:,:,:,1)
            end if
            dt_hydro_block(blk%blk_id)=t_hydro_block(vglobalmax,dx)
            blk=>blk%blk_next
        end do
        time_sys%dt_hydro=minval(dt_hydro_block)
        blk=>blk_head
        do while (associated(blk))
            du=blk%dudt_x*time_sys%dt_hydro
            blk%u=blk%u+du
            call convert_u_to_w_block(blk)
            if (igeometry/=0) then
                call geometry_source(blk,blk%s_geo2)
                call convert_u_to_source_block(blk)
                blk%source=blk%source+(blk%s_geo1+blk%s_geo2)*time_sys%dt_hydro/two
                call convert_source_to_u_block(blk)
                call convert_u_to_w_block(blk)
            end if
            blk=>blk%blk_next
        end do
        nullify(blk)
        deallocate(du,xfluxarray)
    else if (nd==2) then
        !no godunov 2d, use muscl 2d
    end if
    deallocate(dt_hydro_block)
end subroutine godunov_hydro_split

end module godunov
