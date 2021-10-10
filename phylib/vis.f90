module vis
use datastructure
use phylib
use eos
implicit none
contains

subroutine viscous_hydro(blk)
    type(blockdef), pointer :: blk
    if (viscous) then
        call vis_stress_tensor(blk)
        call vis_momentum(blk)
        call vis_energy(blk)
    end if
end subroutine viscous_hydro

subroutine vis_stress_tensor(blk)
    !only need to calculate the viscous tensor on the surfaces of the domain
    !no need to calculate it on the guard cells
    !*******************************careful, SMR need special treatment*********
    type(blockdef), pointer :: blk
    integer :: i,j,ijk(3)
    real(8) :: vrl,vru,vphil,vphiu,rl,ru,r,dr,dphi,nu,rho,mu,vthetal,vthetau,lever,theta,dtheta,thetal, &
        thetau,pos(3),blk_dxyz(3),blk_coords(3),vol,surf_xl,surf_xu,surf_yl,surf_yu
    real(8), dimension(:,:,:), allocatable :: dvrdr,dvrdphi,domegadr,dvphidphi,dvrdtheta,dvthetadtheta, &
        domegadr1,domegadr2,dvphi_sindtheta
    if (igeometry==1) then
        nu=2d-1
        dr=blk%dxyz(1)
        dphi=blk%dxyz(2)
        call allocate_cell_data_block(dvrdr)
        call allocate_cell_data_block(dvrdphi)
        call allocate_cell_data_block(domegadr)
        call allocate_cell_data_block(dvphidphi)
        do j=blk_ylb+1,blk_yub-1
            do i=blk_xlb+1,blk_xub-1
                r=blk%x_center(i)
                dvrdr(i,j,1)=(blk%w(i+1,j,1,2)-blk%w(i-1,j,1,2))/dr/2d0
                dvrdphi(i,j,1)=(blk%w(i,j+1,1,2)-blk%w(i,j-1,1,2))/dphi/2d0
                domegadr(i,j,1)=(blk%w(i+1,j,1,3)/(r+dr)-blk%w(i-1,j,1,3)/(r-dr))/dr/2d0
                dvphidphi(i,j,1)=(blk%w(i,j+1,1,3)-blk%w(i,j-1,1,3))/dphi/2d0
                rl=r-half*dr
                ru=r+half*dr
                vrl=(blk%w(i-1,j,1,2)+blk%w(i,j,1,2))/2d0
                vru=(blk%w(i,j,1,2)+blk%w(i+1,j,1,2))/2d0
                vphil=(blk%w(i,j-1,1,3)+blk%w(i,j,1,3))/2d0
                vphiu=(blk%w(i,j,1,3)+blk%w(i,j+1,1,3))/2d0
                surf_xl=blk%surf1(i-1,j,1)
                surf_xu=blk%surf1(i,j,1)
                surf_yl=blk%surf2(i,j-1,1)
                surf_yu=blk%surf2(i,j,1)
                vol=blk%vol(i,j,1)
                blk%divv(i,j,1)=(surf_xu*vru-surf_xl*vrl+surf_yu*vphiu-surf_yl*vphil)/vol
            end do
        end do
        do j=1,blk_size_ny
            do i=0,blk_size_nx
                rho=min(blk%w(i,j,1,1),blk%w(i+1,j,1,1))
                mu=rho*nu
                r=blk%r_interface(i)
                blk%vis_tensor_xx(i,j,1)=mu*(dvrdr(i,j,1)+dvrdr(i+1,j,1)-(blk%divv(i,j,1)+blk%divv(i+1,j,1))/3d0)
                blk%vis_tensor_xy(i,j,1)=mu*((dvrdphi(i,j,1)/blk%x_center(i)+dvrdphi(i+1,j,1)/blk%x_center(i+1))/2d0    &
                    +(blk%x_center(i)*domegadr(i,j,1)+blk%x_center(i+1)*domegadr(i+1,j,1))/2d0)
            end do
        end do
        do j=0,blk_size_ny
            do i=1,blk_size_nx
                r=blk%x_center(i)
                rho=min(blk%w(i,j,1,1),blk%w(i,j+1,1,1))
                mu=rho*nu
                blk%vis_tensor_yx(i,j,1)=mu*((dvrdphi(i,j,1)+dvrdphi(i,j+1,1))/r/2d0+r*(domegadr(i,j,1)+domegadr(i,j+1,1))/2d0)
                blk%vis_tensor_yy(i,j,1)=mu*((dvphidphi(i,j,1)+dvphidphi(i,j+1,1)+blk%w(i,j,1,2)+blk%w(i,j+1,1,2))/r   &
                    -(blk%divv(i,j,1)+blk%divv(i,j+1,1)/3d0))
            end do
        end do
        deallocate(dvrdr,dvrdphi,domegadr,dvphidphi)
    else if (igeometry==2) then
        nu=1d-1
        blk_dxyz=blk%dxyz
        blk_coords=blk%pos
        dr=blk%dxyz(1)
        dtheta=blk%dxyz(2)
        call allocate_cell_data_block(dvrdr)
        call allocate_cell_data_block(dvrdtheta)
        call allocate_cell_data_block(domegadr1)
        call allocate_cell_data_block(domegadr2)
        call allocate_cell_data_block(dvthetadtheta)
        call allocate_cell_data_block(dvphi_sindtheta)
        do j=blk_ylb+1,blk_yub-1
            do i=blk_xlb+1,blk_xub-1
                ijk=(/i,j,1/)
                call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                r=pos(1)
                theta=pos(2)
                lever=blk%lever(i,j,1)
                dvrdr(i,j,1)=(blk%w(i+1,j,1,2)-blk%w(i-1,j,1,2))/dr/2d0
                dvrdtheta(i,j,1)=(blk%w(i,j+1,1,2)-blk%w(i,j-1,1,2))/dtheta/2d0
                domegadr1(i,j,1)=(blk%w(i+1,j,1,3)/blk%x_center(i+1)-blk%w(i-1,j,1,3)/blk%x_center(i-1))/dr/2d0
                domegadr2(i,j,1)=(blk%vphi(i+1,j,1)/blk%x_center(i+1)-blk%vphi(i-1,j,1)/blk%x_center(i-1))/dr/2d0
                dvthetadtheta(i,j,1)=(blk%w(i,j+1,1,3)-blk%w(i,j-1,1,3))/dtheta/2d0
                dvphi_sindtheta(i,j,1)=(blk%vphi(i,j+1,1)/sin(theta+dtheta)-blk%vphi(i,j-1,1)/sin(theta-dtheta))/dtheta/2d0
                rl=r-half*dr
                ru=r+half*dr
                thetal=theta-half*dtheta
                thetau=theta+half*dtheta
                vrl=(blk%w(i-1,j,1,2)+blk%w(i,j,1,2))/2d0
                vru=(blk%w(i,j,1,2)+blk%w(i+1,j,1,2))/2d0
                vthetal=(blk%w(i,j-1,1,3)+blk%w(i,j,1,3))/2d0
                vthetau=(blk%w(i,j,1,3)+blk%w(i,j+1,1,3))/2d0
                surf_xl=blk%surf1(i-1,j,1)
                surf_xu=blk%surf1(i,j,1)
                surf_yl=blk%surf2(i,j-1,1)
                surf_yu=blk%surf2(i,j,1)
                vol=blk%vol(i,j,1)
                blk%divv(i,j,1)=(vru*surf_xu-vrl*surf_xl+vthetau*surf_yu-vthetal*surf_yl)/vol
            end do
        end do
        do j=1,blk_size_ny
            do i=0,blk_size_nx
                rho=min(blk%w(i,j,1,1),blk%w(i+1,j,1,1))
                mu=rho*nu
                blk%vis_tensor_xx(i,j,1)=mu*(dvrdr(i,j,1)+dvrdr(i+1,j,1)-(blk%divv(i,j,1)+blk%divv(i+1,j,1))/3d0)
                blk%vis_tensor_xy(i,j,1)=mu*(dvrdtheta(i,j,1)/blk%x_center(i)+dvrdtheta(i+1,j,1)/blk%x_center(i+1)  &
                    +blk%x_center(i)*domegadr1(i,j,1)+blk%x_center(i+1)*domegadr1(i+1,j,1))/2d0
                blk%vis_tensor_rphi(i,j,1)=mu*(blk%x_center(i)*domegadr2(i,j,1)+blk%x_center(i+1)*domegadr2(i+1,j,1))/2d0
            end do
        end do
        do j=0,blk_size_ny
            do i=1,blk_size_nx
                r=blk%x_center(i)
                rho=min(blk%w(i,j,1,1),blk%w(i,j+1,1,1))
                mu=rho*nu
                blk%vis_tensor_yx(i,j,1)=mu*(dvrdtheta(i,j,1)/r+dvrdtheta(i,j+1,1)/r+r*domegadr1(i,j,1)+r*domegadr1(i,j+1,1))/2d0
                blk%vis_tensor_yy(i,j,1)=mu*(dvthetadtheta(i,j,1)/r+dvthetadtheta(i,j+1,1)/r    &
                    +blk%w(i,j,1,2)/r+blk%w(i,j+1,1,2)/r-(blk%divv(i,j,1)+blk%divv(i,j+1,1))/3d0)
                blk%vis_tensor_thetaphi(i,j,1)=mu*(blk%lever(i,j,1)/r**2*dvphi_sindtheta(i,j,1) &
                    +blk%lever(i,j+1,1)/r**2*dvphi_sindtheta(i,j+1,1))/2d0
            end do
        end do
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                ijk=(/i,j,1/)
                call ijk_to_coords(ijk,blk_dxyz,blk_coords,pos)
                r=pos(1)
                theta=pos(2)
                rho=blk%w(i,j,1,1)
                mu=rho*nu
                blk%vis_tensor_phiphi(i,j,1)=mu*2d0*(blk%w(i,j,1,2)/blk%x_center(i)+blk%w(i,j,1,3)*cos(theta)/blk%lever(i,j,1)-blk%divv(i,j,1)/3d0)
                blk%vis_tensor_phitheta(i,j,1)=mu*blk%lever(i,j,1)/r**2*dvphi_sindtheta(i,j,1)
            end do
        end do
        deallocate(dvrdr,dvrdtheta,domegadr1,domegadr2,dvthetadtheta,dvphi_sindtheta)
    end if
end subroutine vis_stress_tensor

subroutine vis_momentum(blk)
    !only need to compute the viscous force within the domain cell
    type(blockdef), pointer :: blk
    integer :: i,j
    real(8) :: txxl,txxu,txyl,txyu,tyxl,tyxu,tyyl,tyyu,rl,ru,r,dr,dt,rho,dphi,vol,surf_xl,surf_xu,surf_yl,surf_yu,  &
        am1,am2,tzz,tzy,surf_z,txzl,txzu,tyzl,tyzu
    dt=time_sys%dt_hydro
    if (igeometry==1) then
        dr=blk%dxyz(1)
        dphi=blk%dxyz(2)
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                vol=blk%vol(i,j,1)
                surf_xl=blk%surf1(i-1,j,1)
                surf_xu=blk%surf1(i,j,1)
                surf_yl=blk%surf2(i,j-1,1)
                surf_yu=blk%surf2(i,j,1)
                txxl=blk%vis_tensor_xx(i-1,j,1)
                txxu=blk%vis_tensor_xx(i,j,1)
                txyl=blk%vis_tensor_xy(i-1,j,1)
                txyu=blk%vis_tensor_xy(i,j,1)
                tyxl=blk%vis_tensor_yx(i,j-1,1)
                tyxu=blk%vis_tensor_yx(i,j,1)
                tyyl=blk%vis_tensor_yy(i,j-1,1)
                tyyu=blk%vis_tensor_yy(i,j,1)
                r=blk%x_center(i)
                rl=r-half*dr
                ru=r+half*dr
                blk%w(i,j,1,2)=(blk%u(i,j,1,2)*vol+((txxu*surf_xu-txxl*surf_xl)+(tyxu*surf_yu-tyxl*surf_yl)     &
                    -(tyyl+tyyu)/2*(surf_xu-surf_xl))*dt)/vol/blk%w(i,j,1,1)
                am1=blk%u(i,j,1,3)*r
                am2=(am1*vol+((ru*surf_xu*txyu-rl*surf_xl*txyl)+(r*surf_yu*tyyu-r*surf_yl*tyyl))*dt)/vol
                blk%w(i,j,1,3)=am2/r/blk%w(i,j,1,1)
            end do
        end do
        call convert_w_to_u_block(blk)
    else if (igeometry==2) then
        do j=1,blk_size_ny
            do i=1,blk_size_nx
                r=blk%x_center(i)
                vol=blk%vol(i,j,1)
                surf_xl=blk%surf1(i-1,j,1)
                surf_xu=blk%surf1(i,j,1)
                surf_yl=blk%surf2(i,j-1,1)
                surf_yu=blk%surf2(i,j,1)
                surf_z=blk%surf_phi(i,j,1)
                txxl=blk%vis_tensor_xx(i-1,j,1)
                txxu=blk%vis_tensor_xx(i,j,1)
                txyl=blk%vis_tensor_xy(i-1,j,1)
                txyu=blk%vis_tensor_xy(i,j,1)
                txzl=blk%vis_tensor_rphi(i-1,j,1)
                txzu=blk%vis_tensor_rphi(i,j,1)
                tyxl=blk%vis_tensor_yx(i,j-1,1)
                tyxu=blk%vis_tensor_yx(i,j,1)
                tyyl=blk%vis_tensor_yy(i,j-1,1)
                tyyu=blk%vis_tensor_yy(i,j,1)
                tyzl=blk%vis_tensor_thetaphi(i,j-1,1)
                tyzu=blk%vis_tensor_thetaphi(i,j,1)
                tzz=blk%vis_tensor_phiphi(i,j,1)
                tzy=blk%vis_tensor_phitheta(i,j,1)
                !blk%w(i,j,1,2)=(blk%u(i,j,1,2)*vol+((txxu*surf_xu-txxl*surf_xl)+(tyxu*surf_yu-tyxl*surf_yl)&
                !    -((tyyl+tyyu)/2+tzz)*(surf_xu-surf_xl)/2)*dt)/vol/blk%w(i,j,1,1)
                !blk%w(i,j,1,3)=(blk%u(i,j,1,3)*r*vol+((ru*txyu*surf_xu-rl*txyl*surf_xl)+r*(tyyu*surf_yu-tyyl*surf_yl)&
                !    -r*tzz*(surf_yu-surf_yl))*dt)/vol/r/blk%w(i,j,1,1)
                am1=blk%w(i,j,1,1)*blk%vphi(i,j,1)*blk%lever(i,j,1)
                am2=(am1*vol+(blk%lever_r(i,j,1)*txzu*surf_xu-blk%lever_r(i-1,j,1)*txzl*surf_xl &
                    +blk%lever_theta(i,j,1)*tyzu*surf_yu-blk%lever_theta(i,j-1,1)*tyzl*surf_yl)*dt)/vol
                blk%vphi(i,j,1)=am2/blk%w(i,j,1,1)/blk%lever(i,j,1)
            end do
        end do
        call convert_w_to_u_block(blk)
    end if
end subroutine vis_momentum

subroutine vis_energy(blk)
    !viscouse heating
    type(blockdef), pointer :: blk
    if (igeometry==1) then
    else if (igeometry==2) then
    end if
end subroutine vis_energy

end module vis
