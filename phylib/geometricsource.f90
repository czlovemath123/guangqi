module geometricsource
use datastructure
use phylib
use eos
implicit none
contains

subroutine geometry_source(blk,s_geo)
    !conservative source
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:,:), allocatable :: s_geo
    if (nd==1) then
        if (igeometry==1) then
            !call cylindrical1d_block(blk,s_geo)
        else if (igeometry==2) then
            call polar1d_block(blk,s_geo)
        end if
    else if (nd==2) then
    end if
end subroutine geometry_source

subroutine polar1d_block(blk,s_geo)
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:,:), allocatable :: s_geo
    integer :: i,ijk(3)
    real(8) :: p
    do i=blk_xlb,blk_xub
        p=blk%w(i,1,1,5)
        s_geo(i,1,1,1:5)=0d0
        s_geo(i,1,1,2)=p*(blk%surf1(i,1,1)-blk%surf1(i-1,1,1))
    end do
end subroutine polar1d_block

subroutine polar2d_geometry_predict(blk)
    !called in the corrector step of 2d spherical coordinate, include gravity
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:,:), allocatable :: s_geo
    integer :: i,j,ijk(3)
    real(8) :: p,theta,r,dxyz(3),rho,vol,vr,vtheta,vphi,g
    dxyz=blk%dxyz
    do j=blk_ylb+1,blk_yub-1
        do i=blk_xlb+1,blk_xub-1
            vol=blk%vol(i,j,1)
            rho=blk%predict_w(i,j,1,1)
            vr=blk%predict_w(i,j,1,2)
            vtheta=blk%predict_w(i,j,1,3)
            vphi=blk%vphi(i,j,1)
            !r=blk%pos(1)+dxyz(1)*(i-half)
            r=blk%x_center(i)
            g=gr*central_star%core%mass/(blk%pos(1)+dxyz(1)*(i-half))**2
            blk%s_geo1(i,j,1,1:5)=0d0
            blk%s_geo1(i,j,1,2)=rho*(vtheta**2+vphi**2)/r-rho*g
            blk%s_geo1(i,j,1,3)=rho*vphi**2*(blk%surf2(i,j,1)-blk%surf2(i,j-1,1))/vol-rho*vr*vtheta/r
        end do
    end do
end subroutine polar2d_geometry_predict

subroutine polar2d_geometry_correct(blk)
    !called in the corrector step of 2d spherical coordinate, include gravity
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:,:), allocatable :: s_geo
    integer :: i,j,ijk(3)
    real(8) :: p,ctheta,r,dxyz(3),rho,vol,vr,vtheta,vphi,g,val,vphil,vphir,flux1,flux2,r1,r2
    dxyz=blk%dxyz
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            vol=blk%vol(i,j,1)
            p=blk%w(i,j,1,5)
            rho=blk%w(i,j,1,1)
            vr=blk%w(i,j,1,2)
            vtheta=blk%w(i,j,1,3)
            vphi=blk%vphi(i,j,1)
            !r=blk%pos(1)+dxyz(1)*(i-half)
            r=blk%x_center(i)
            flux1=blk%xflux(i-1,j,1,3)
            flux2=blk%xflux(i,j,1,3)
            r1=blk%r_interface(i-1)
            r2=blk%r_interface(i)
            g=gr*central_star%core%mass/r
            !g=gr*central_star%core%mass*3d0/(r1**2+r1*r2+r2**2)
            !g=gr*central_star%core%mass/(blk%pos(1)+dxyz(1)*(i-half))**2d0
            !vphi=0d0;g=0d0
            !vtheta=0d0;vphi=0d0
            blk%s_geo2(i,j,1,1:5)=0d0
            blk%s_geo2(i,j,1,2)=0.5d0*rho*(vtheta**2d0+vphi**2d0)*(blk%surf1(i,j,1)-blk%surf1(i-1,j,1))/vol-0.5*rho*g*(blk%surf1(i,j,1)-blk%surf1(i-1,j,1))/vol
            !print *,'2',blk%s_geo2(i,j,1,2),0.5d0*rho*(vtheta**2d0+vphi**2d0)*(blk%surf1(i,j,1)-blk%surf1(i-1,j,1))/vol,rho*g
            blk%s_geo2(i,j,1,3)=rho*vphi**2d0*(blk%surf2(i,j,1)-blk%surf2(i,j-1,1))/vol &
                -0.25d0*(flux2+flux1)*(blk%surf1(i,j,1)-blk%surf1(i-1,j,1))/vol
            !print *,'3',rho*vphi**2d0*(blk%surf2(i,j,1)-blk%surf2(i,j-1,1))/vol,0.25d0*(flux2+flux1)*(blk%surf1(i,j,1)-blk%surf1(i-1,j,1))/vol
            !print *,blk%s_geo2(i,j,1,1:5)
            !print *,blk%s_geo2(i,j,1,2)/blk%w(i,j,1,1),vtheta,vphi**2d0/r,g
            !print *,blk%s_geo2(i,j,1,3)/blk%w(i,j,1,1)
            !stop
            !print *,rho*((vtheta**2d0+vphi**2d0)/r-g)+(p*blk%surf1(i,j,1)-p*blk%surf1(i-1,j,1))/vol
            !if (blk%loc_type==9.and.i==1.and.j==1) then
            !    print *,'9',r,p,vol,blk%s_geo2(i,j,1,3)
            !else if (blk%loc_type==3.and.i==1.and.j==16) then
            !    print *,'3',r,p,vol,blk%s_geo2(i,j,1,3)
            !end if
            !val=(rho*(vtheta**2+vphi**2)/r-rho*g)*vol/(p*(blk%surf1(i,j,1)-blk%surf1(i-1,j,1)))
            !if (val>1d0) print *,val,blk%loc_type,i,j
            !print *,blk%s_geo2(i,j,1,2)*time_sys%dt_hydro/vol/blk%w(i,j,1,1),   &
            !    p*(blk%surf1(i,j,1)-blk%surf1(i-1,j,1))*time_sys%dt_hydro/vol/blk%w(i,j,1,1),   &
            !    (rho*(vtheta**2+vphi**2)/r-rho*g)*time_sys%dt_hydro/blk%w(i,j,1,1)
        end do
    end do
    !stop
end subroutine polar2d_geometry_correct

!subroutine cylindrical2d_geometry_predict(blk)
!    !called in the corrector step of 2d cylindrical coordinate, include gravity
!    type(blockdef), pointer :: blk
!    integer :: i,j
!    real(8) :: p,rho,vr,vphi,dxyz(3),r,g,vol
!    dxyz=blk%dxyz
!    do j=blk_ylb+1,blk_yub-1
!        do i=blk_xlb+1,blk_xub-1
!            vol=blk%vol(i,j,1)
!            p=blk%predict_w(i,j,1,5)
!            rho=blk%predict_w(i,j,1,1)
!            vr=blk%predict_w(i,j,1,2)
!            vphi=blk%predict_w(i,j,1,3)
!            r=blk%pos(1)+dxyz(1)*(i-half)
!            g=gr*central_star%core%mass/r**2
!            blk%s_geo1(i,j,1,1:5)=0d0
!            blk%s_geo1(i,j,1,2)=p*(blk%surf1(i,j,1)-blk%surf1(i-1,j,1))+rho*(vphi**2/r-g)*vol
!        end do
!    end do
!end subroutine cylindrical2d_geometry_predict

subroutine cylindrical2d_geometry_correct(blk)
    !called in the corrector step of 2d cylindrical coordinate, include gravity
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: p,rho,vr,vphi,dxyz(3),r,g,vol
    dxyz=blk%dxyz
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            vol=blk%vol(i,j,1)
            p=blk%w(i,j,1,5)
            rho=blk%w(i,j,1,1)
            vr=blk%w(i,j,1,2)
            vphi=blk%w(i,j,1,3)
            r=blk%pos(1)+dxyz(1)*(i-half)
            g=gr*central_star%core%mass/r**2
            blk%s_geo2(i,j,1,1:5)=0d0
            blk%s_geo2(i,j,1,2)=p*(blk%surf1(i,j,1)-blk%surf1(i-1,j,1))+rho*(vphi**2/r-g)*vol
        end do
    end do
end subroutine cylindrical2d_geometry_correct

end module geometricsource
