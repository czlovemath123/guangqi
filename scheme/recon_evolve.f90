module recon_evolve
use mathlib
use datastructure
use eos
use phylib
use limiters
use source_control

implicit none

contains

subroutine reconstruct_hydro(blk)
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:,:), allocatable :: w,wl,wr
    real(8), dimension(:,:,:), allocatable :: temp
    real(8), dimension(5) :: v1,v2,v3,slp_left,slp_right,slp_central,slp
    real(8) :: xl,xc,xr,xl_face,xr_face,yl,yc,yr,yl_face,yr_face
    integer :: i,j,k
    if (nd==1) then
        do i=blk_xlb+1,blk_xub-1
            v1=blk%w(i-1,1,1,1:5)
            v2=blk%w(i,1,1,1:5)
            v3=blk%w(i+1,1,1,1:5)
            xl=blk%x_center(i-1)
            xc=blk%x_center(i)
            xr=blk%x_center(i+1)
            slp_left=(v2-v1)/(xc-xl)
            slp_right=(v3-v2)/(xr-xc)
            slp_central=(v3-v1)/(xr-xl)
            do k=1,5
                slp(k)=find_the_slope(slp_left(k),slp_right(k),slp_central(k))
            end do
            blk%xslp(i,1,1,1:5)=slp
        end do
    else if (nd==2) then
        do j=1,blk_size_ny
            do i=blk_xlb+1,blk_xub-1
                v1=blk%w(i-1,j,1,1:5)
                v2=blk%w(i,j,1,1:5)
                v3=blk%w(i+1,j,1,1:5)
                xl=blk%x_center(i-1)
                xc=blk%x_center(i)
                xr=blk%x_center(i+1)
                slp_left=(v2-v1)/(xc-xl)
                slp_right=(v3-v2)/(xr-xc)
                slp_central=(v3-v1)/(xr-xl)
                do k=1,5
                    slp(k)=find_the_slope(slp_left(k),slp_right(k),slp_central(k))
                end do
                xl_face=blk%x_interface(i-1)
                xr_face=blk%x_interface(i)
                blk%xslp(i,j,1,1:5)=slp
                blk%w_xl(i,j,1,1:5)=blk%w(i,j,1,1:5)+(xl_face-xc)*slp
                blk%w_xr(i,j,1,1:5)=blk%w(i,j,1,1:5)+(xr_face-xc)*slp
            end do
        end do
        do j=blk_ylb+1,blk_yub-1
            do i=1,blk_size_nx
                v1=blk%w(i,j-1,1,1:5)
                v2=blk%w(i,j,1,1:5)
                v3=blk%w(i,j+1,1,1:5)
                yl=blk%y_center(j-1)
                yc=blk%y_center(j)
                yr=blk%y_center(j+1)
                slp_left=(v2-v1)/(yc-yl)
                slp_right=(v3-v2)/(yr-yc)
                slp_central=(v3-v1)/(yr-yl)
                do k=1,5
                    slp(k)=find_the_slope(slp_left(k),slp_right(k),slp_central(k))
                end do
                yl_face=blk%y_interface(j-1)
                yr_face=blk%y_interface(j)
                blk%yslp(i,j,1,1:5)=slp
                blk%w_yl(i,j,1,1:5)=blk%w(i,j,1,1:5)+(yl_face-yc)*slp
                blk%w_yr(i,j,1,1:5)=blk%w(i,j,1,1:5)+(yr_face-yc)*slp
            end do
        end do
    end if
end subroutine reconstruct_hydro

subroutine evolve_block(blk)
    !used only in 1d
    type(blockdef), pointer :: blk
    real(8) :: dwl(5),dwr(5),w(5),dw(5),wl(5),wr(5),slp(5),egv,temp,cs,xl,xc,xr,dt
    integer :: i,j,k,ijk(3)
    dt=time_sys%dt_hydro
    do i=blk_xlb+1,blk_xub-1
        dwl=0d0;dwr=0d0;dw=0d0
        w=blk%w(i,1,1,1:5)
        xl=blk%x_interface(i-1)
        xc=blk%x_center(i)
        xr=blk%x_interface(i)
#if     ieos==1
        cs=sqrt(gamma_gas*w(5)/w(1))
#elif   ieos==2
        temp=blk%temp(i,1,1)
        cs=calculate_adiabatic_cs(w(1),temp)
#endif
        slp=blk%xslp(i,1,1,1:5)
        dwl(1)=(xl-xc)*slp(1)-dt/2*(w(2)*slp(1)+w(1)*slp(2))
        dwl(2)=(xl-xc)*slp(2)-dt/2*(w(2)*slp(2)+slp(5)/w(1))
        dwl(3)=0d0
        dwl(4)=0d0
        dwl(5)=(xl-xc)*slp(5)-dt/2*(w(1)*cs**2*slp(2)+w(2)*slp(5))
        dwr(1)=(xr-xc)*slp(1)-dt/2*(w(2)*slp(1)+w(1)*slp(2))
        dwr(2)=(xr-xc)*slp(2)-dt/2*(w(2)*slp(2)+slp(5)/w(1))
        dwr(3)=0d0
        dwr(4)=0d0
        dwr(5)=(xr-xc)*slp(5)-dt/2*(w(1)*cs**2*slp(2)+w(2)*slp(5))
        if (igeometry==2) then
            egv=blk%egv(i,1,1)
            call evolve_geometric_source(w,egv,xc,dw)
            dwl=dwl+dw*dt/2
            dwr=dwr+dw*dt/2
        end if
        blk%w_xl(i,1,1,1:5)=w+dwl
        blk%w_xr(i,1,1,1:5)=w+dwr
#if     ieos==2
        wl=w+dwl
        wr=w+dwr
        blk%temp_xl(i,1,1)=calculate_temp_from_p_rho(wl(5),wl(1))
        blk%temp_xr(i,1,1)=calculate_temp_from_p_rho(wr(5),wr(1))
        blk%egv_xl(i,1,1)=egvrhot(wl(1),blk%temp_xl(i,1,1))
        blk%egv_xr(i,1,1)=egvrhot(wr(1),blk%temp_xr(i,1,1))
#endif
    end do
end subroutine evolve_block

subroutine evolve_geometric_source(w,egv,r,dw)
    real(8) :: w(5),egv,dw(5),r
    dw=0d0
    dw(1)=-2d0/r*w(1)*w(2)
    dw(5)=-2d0/r*(w(5)+egv)*w(2)
end subroutine evolve_geometric_source

end module recon_evolve
