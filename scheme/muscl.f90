module muscl
use mathlib
use phylib
use hydro
use eos
use datastructure
use communication
!use limiters
use recon_evolve
use boundary
implicit none

contains

subroutine vanleer_hydro_unsplit()
    !based on Stone & Gardiner 2009 "A simple unsplit Godunov method for multidimensional MHD"
    if (nd==1) then
        call collective_sub(estimate_block_dt_hydro,collective_processor_min,collective_rank0_min,time_sys%dt_hydro)
        call muscl_hydro_split()
    else if (nd==2) then
        call collective_sub(estimate_block_dt_hydro,collective_processor_min,collective_rank0_min,time_sys%dt_hydro)
        call blk_traversal(vanleer_predictor)
        !call check_symmetry()
        call communicate_hydro()
        call applyboundconds()
        call blk_traversal(vanleer_corrector)
    end if
end subroutine vanleer_hydro_unsplit

subroutine check_symmetry()
    type(blockdef), pointer :: blk1,blk2
    real(8) :: c1,c2
    integer :: i,j
    blk1=>blk_processor_head
    blk2=>blk_processor_head%blk_next
    print *,'ntimestep=',time_sys%ntimestep
    do j=1,2
        if (j<=2) then
            print *,'blk1 id=',blk1%blk_id,'blk2 id=',blk2%blk_id
            !print *,norm2(blk1%xflux(blk_size_nx,1:blk_size_ny,1,1:5)-blk2%xflux(0,1:blk_size_ny,1,1:5))
            !print *,blk1%xflux(blk_size_nx-1,1:blk_size_ny,1,1:5)-blk2%xflux(1,1:blk_size_ny,1,1:5)
            !print *,norm2(blk1%w(blk_size_nx,1:blk_size_ny,1,1:5)-blk2%w(1,1:blk_size_ny,1,1:5))
            c1=norm2(blk1%w(blk_size_nx,1:blk_size_ny,1,1)-blk2%w(1,1:blk_size_ny,1,1))
            c2=norm2(blk1%w(blk_size_nx,1:blk_size_ny,1,1)-blk2%w(1,1:blk_size_ny,1,1))
            if (c1>=epsilon(1d0).or.c2>=epsilon(1d0)) then
                print *,c1,c2
                !print *,blk1%w(blk_size_nx,1:blk_size_ny,1,1)-blk2%w(1,1:blk_size_ny,1,1)
                !print *,blk1%w(blk_size_nx,1:blk_size_ny,1,5)
                !print *,blk2%w(1,1:blk_size_ny,1,5)
                !stop
            end if
            !do i=1,blk_size_ny
            !    !print *,i,blk1%w(blk_size_nx-4,i,1,1:5)-blk2%w(5,i,1,1:5),blk1%w(blk_size_nx-4,i,1,2)+blk2%w(5,i,1,2)
            !    !print *,blk1%xflux(blk_size_nx-1,i,1,2)-blk2%xflux(1,i,1,2)
            !    !print *,blk1%yflux(blk_size_nx,i,1,2)+blk2%yflux(1,i,1,2)
            !    !print *,blk1%yflux(blk_size_nx,i-1,1,2)+blk2%yflux(1,i-1,1,2)
            !    !print *,blk1%xflux(blk_size_nx-5,i,1,1:5)-blk1%xflux(blk_size_nx-4,i,1,1:5)+blk1%yflux(blk_size_nx-4,i-1,1,1:5)-blk1%yflux(blk_size_nx-4,i,1,1:5)
            !    !print *,blk2%xflux(4,i,1,1:5)-blk2%xflux(5,i,1,1:5)+blk2%yflux(5,i-1,1,1:5)-blk2%yflux(5,i,1,1:5)
            !    print *,i
            !    print *,blk1%xflux(blk_size_nx,i,1,1:5)
            !    !print *,blk1%xflux(blk_size_nx-4,i,1,2)
            !    !print *,blk1%yflux(blk_size_nx-4,i-1,1,2)-blk1%yflux(blk_size_nx-4,i,1,2)
            !    print *,blk2%xflux(0,i,1,1:5)
            !    !print *,blk2%xflux(5,i,1,2)
            !    !print *,blk2%yflux(5,i-1,1,2)-blk2%yflux(5,i,1,2)
            !    !print *,blk1%xflux(blk_size_nx-4,i,1,1:5)+blk2%xflux(4,i,1,1:5)
            !    !print *,blk1%xflux(blk_size_nx-5,i,1,1:5)+blk2%xflux(5,i,1,1:5)
            !    !print *,blk1%w_xl(blk_size_nx-4,i,1,1:5)
            !    !print *,blk2%w_xr(5,i,1,1:5)
            !    !print *,blk1%w_xr(blk_size_nx-5,i,1,1:5)
            !    !print *,blk2%w_xl(6,i,1,1:5)
            !    !print *,blk1%w0(blk_size_nx-3,i,1,1:5)
            !    !print *,blk2%w0(4,i,1,1:5)
            !    !print *,blk1%w0(blk_size_nx-2,i,1,1:5)
            !    !print *,blk2%w0(3,i,1,1:5)
            !    !print *,blk1%w0(blk_size_nx-4,i,1,1:5)
            !    !print *,blk2%w0(5,i,1,1:5)
            !end do
        end if
        if (j/=2) then
            blk1=>blk2%blk_next
            blk2=>blk1%blk_next
        end if
    end do
end subroutine check_symmetry

subroutine estimate_block_dt_hydro(blk,blk_result)
    !estimate dt_hydro on all interfaces, save the minimum value to blk_result
    !t=0 timestep dt estimate is different
    type(blockdef), pointer :: blk
    real(8), dimension(:,:,:), allocatable :: varray,varray_x,varray_y
    real(8) :: vmax,cs,rho,temp,dt(nd),vblockmax(nd)
    real(8), allocatable :: dty(:)
    character(len=1) :: dir
    integer :: i,j,k,ierr
    real(8) :: blk_result
    if (nd==1) then
        if (time_sys%ntimestep==0.and.time_sys%t==0) then
            blk%w_xl=blk%w;blk%w_xr=blk%w
#if       ieos==2
            blk%temp_xl=blk%temp;blk%temp_xr=blk%temp
            blk%egv_xl=blk%egv;blk%egv_xr=blk%egv
#endif
            call hllc_muscl(blk)
            blk_result=blk%dt_hydro_blk
        else
            call allocate_cell_data_block(varray)
            do i=blk_xlb,blk_xub
                rho=blk%w(i,1,1,1)
                temp=blk%temp(i,1,1)
                cs=calculate_adiabatic_cs(rho,temp)
                varray(i,1,1)=abs(blk%w(i,1,1,2))+cs
            end do
            vmax=maxval(varray)
            blk_result=blk%dxyz(1)*CFL/vmax
            deallocate(varray)
        end if
    else if (nd==2) then
        if (time_sys%ntimestep==0.and.time_sys%t==0) then
            blk%w_xl=blk%w;blk%w_xr=blk%w;blk%w_yl=blk%w;blk%w_yr=blk%w
#if       ieos==2
            blk%temp_xl=blk%temp;blk%temp_xr=blk%temp;blk%temp_yl=blk%temp;blk%temp_yr=blk%temp
            blk%egv_xl=blk%egv;blk%egv_xr=blk%egv;blk%egv_yl=blk%egv;blk%egv_yr=blk%egv
#endif
            call hllc_muscl(blk)
            blk%w_xl=0d0;blk%w_xr=0d0;blk%w_yl=0d0;blk%w_yr=0d0
            blk%xflux=0d0;blk%yflux=0d0
            blk_result=blk%dt_hydro_blk
        else
            call allocate_cell_data_block(varray_x)
            call allocate_cell_data_block(varray_y)
            do j=blk_ylb,blk_yub
                do i=blk_xlb,blk_xub
                    rho=blk%w(i,j,1,1)
                    temp=blk%temp(i,j,1)
                    cs=calculate_adiabatic_cs(rho,temp)
                    varray_x(i,j,1)=abs(blk%w(i,j,1,2))+cs
                    varray_y(i,j,1)=abs(blk%w(i,j,1,3))+cs
                end do
            end do
            if (igeometry==0) then
                vmax=max(maxval(varray_x),maxval(varray_y))
                blk_result=blk%dxyz(1)*CFL/vmax
            else if (igeometry==1.or.igeometry==2) then
                dt(1)=blk%dxyz(1)*CFL/maxval(varray_x)
                allocate(dty(blk_size_nx))
                do i=1,blk_size_nx
                    dty(i)=blk%dy(i)*CFL/maxval(varray_y(i,:,1))
                end do
                blk_result=min(dt(1),minval(dty))
                deallocate(dty)
            end if
            deallocate(varray_x,varray_y)
        end if
    end if
end subroutine estimate_block_dt_hydro

subroutine muscl_hydro_split()
    type(blockdef), pointer :: blk
    real(8) :: wl(5),wr(5),xflux(5),dx,egv,u(5),w(5),source(5),dvr1,dvr2,mflux_l,mflux_r,s1,s2,vol,     &
        cs,mach,rho,temp,m,dxyz(3),ek,phil,phir,phi,vr_speed,vr,vr_estimate,dt,ek_old,   &
        ek_rank,ek_total,gp_rank,gp_total,ek_pre,gp_pre,e_in,eg_rank,eg_total,eg_pre,p
    integer :: i,j,level,ierr
    logical :: energy_eq
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        call reconstruct_hydro(blk)
        call evolve_block(blk)
        call hllc_muscl(blk)
        if (iradiation/=0) then
            blk%mass_flux_x=0d0!blk%xflux(:,:,:,1)
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    call communicate_flux()
    dt=time_sys%dt_hydro
    blk=>blk_processor_head
    do i=1,np_nblk(rank+1)
        !hydrodynamical conservation law
        do j=1,blk_size_nx
            s1=blk%surf1(j-1,1,1)
            s2=blk%surf1(j,1,1)
            vol=blk%vol(j,1,1)
            blk%u(j,1,1,1:5)=blk%u(j,1,1,1:5)+(blk%xflux(j-1,1,1,1:5)*s1-blk%xflux(j,1,1,1:5)*s2)*dt/vol
        end do
        if (igeometry==2.and.igravity==1) then
            !use energy equation to recover the momentum equation
            m=central_star%core%mass
            dx=blk%dxyz(1)
            do j=1,blk_size_nx
                vol=blk%vol(j,1,1)
                !use the previous step mach number to determine the gravity treatment
                rho=blk%w(j,1,1,1)
                temp=blk%temp(j,1,1)
                cs=calculate_adiabatic_cs(rho,temp)
                mach=blk%w(j,1,1,2)/cs
                if (mach>5.and.abs(blk%w(j,1,1,2))>1d6.and.energy_conserve_formalism) then
                    u=blk%u(j,1,1,1:5)
                    call convert_u_to_w_cell(u,w,temp,egv)
                    ek=u(5)-egv
                    mflux_l=blk%xflux(j-1,1,1,1)
                    mflux_r=blk%xflux(j,1,1,1)
                    if (mflux_l>0) then
                        if (j==1) then
                            if (associated(blk,blk_head)) then
                                phil=blk%gpotential(0,1,1)
                            else
                                if (allocated(blk%blk_pre%gpotential)) then
                                    phil=blk%blk_pre%gpotential(blk_size_nx,1,1)
                                else
                                    phil=-m*gr/(blk%blk_pre%pos(1)+(blk_size_nx-half)*blk%blk_pre%dxyz(1))
                                end if
                            end if
                        else
                            phil=blk%gpotential(j-1,1,1)
                        end if
                    else
                        phil=blk%gpotential(j,1,1)
                    end if
                    if (mflux_r>0) then
                        phir=blk%gpotential(j,1,1)
                    else
                        if (j==blk_size_nx) then
                            if (associated(blk,blk_tail)) then
                                phir=blk%gpotential(blk_size_nx+1,1,1)
                            else
                                if (allocated(blk%blk_next%gpotential)) then
                                    phir=blk%blk_next%gpotential(1,1,1)
                                else
                                    phir=-m*gr/(blk%blk_next%pos(1)+half*blk%blk_next%dxyz(1))
                                end if
                            end if
                        else
                            phir=blk%gpotential(j+1,1,1)
                        end if
                    end if
                    phi=blk%gpotential(j,1,1)
                    s1=blk%surf1(j-1,1,1)
                    s2=blk%surf1(j,1,1)
                    ek=ek+(phi*(mflux_r*s2-mflux_l*s1)-(mflux_r*phir*s2-mflux_l*phil*s1))*dt/vol
                    vr_speed=sqrt(2d0*ek/w(1))
                    !u=blk%u(j,1,1,1:5)
                    u(2)=u(2)-4d0*pi*blk%w(j,1,1,1)*gr*m*dx/vol*dt
                    vr_estimate=u(2)/u(1)
                    dvr1=abs(vr_estimate-vr_speed)
                    dvr2=abs(vr_estimate+vr_speed)
                    if (dvr1<dvr2) then
                        vr=vr_speed
                    else
                        vr=-vr_speed
                    end if
                    u(2)=u(1)*vr
                    !u(5)=blk%egv(j,1,1)+ek
                    u(5)=egv+ek
                    call convert_u_to_w_cell(u,w,temp,egv)
                    blk%u(j,1,1,1:5)=u
                    blk%w(j,1,1,1:5)=w
                    blk%temp(j,1,1)=temp
                    blk%egv(j,1,1)=egv
                else
                    u=blk%u(j,1,1,1:5)
                    p=blk%w(j,1,1,5)
                    u(2)=u(2)+p*(blk%surf1(j,1,1)-blk%surf1(j-1,1,1))*dt/vol
                    call convert_u_to_w_cell(u,w,temp,egv)
                    blk%w(j,1,1,1:5)=w
                    blk%temp(j,1,1)=temp
                    blk%egv(j,1,1)=egv
                    call convert_u_to_source_cell(u,source)
                    source(2)=source(2)-4d0*pi*w(1)*gr*m*dx*dt/vol
                    call convert_source_to_u_cell(source,u)
                    blk%u(j,1,1,1:5)=u
                    blk%w(j,1,1,2)=u(2)/u(1)
                end if
            end do
        else if (igeometry==0) then
            call convert_u_to_w_block(blk)
        end if
        if (i/=np_nblk(rank+1)) then
            blk=>blk%blk_next
        end if
    end do
    nullify(blk)
end subroutine muscl_hydro_split

subroutine vanleer_predictor(blk)
    !evolve the state for half dt
    type(blockdef), pointer :: blk
    real(8) :: u(5),xflux1(5),xflux2(5),yflux1(5),yflux2(5),sx1,sx2,sy1,sy2,vol,dt
    integer :: i,j,k
    blk%w0=blk%w;blk%u0=blk%u;blk%temp0=blk%temp;blk%egv0=blk%egv
    blk%w_xl=blk%w;blk%w_xr=blk%w;blk%w_yl=blk%w;blk%w_yr=blk%w
    !call reconstruct_hydro(blk)
    call hllc_muscl(blk)
    dt=time_sys%dt_hydro
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            sx1=blk%surf1(i-1,j,1)
            sx2=blk%surf1(i,j,1)
            sy1=blk%surf2(i,j-1,1)
            sy2=blk%surf2(i,j,1)
            vol=blk%vol(i,j,1)
            u=blk%u0(i,j,1,:)
            xflux1=blk%xflux(i-1,j,1,1:5)
            xflux2=blk%xflux(i,j,1,1:5)
            yflux1=blk%yflux(i,j-1,1,1:5)
            yflux2=blk%yflux(i,j,1,1:5)
            blk%u(i,j,1,1:5)=u+(xflux1*sx1-xflux2*sx2+yflux1*sy1-yflux2*sy2)*half*dt/vol
        end do
    end do
    call convert_u_to_w_block(blk)
end subroutine vanleer_predictor

subroutine vanleer_corrector(blk)
    type(blockdef), pointer :: blk
    real(8) :: u(5),xflux1(5),xflux2(5),yflux1(5),yflux2(5),sx1,sx2,sy1,sy2,vol,dt
    integer :: i,j,k
    call reconstruct_hydro(blk)
    call hllc_muscl(blk)
    dt=time_sys%dt_hydro
    do j=1,blk_size_ny
        do i=1,blk_size_nx
            sx1=blk%surf1(i-1,j,1)
            sx2=blk%surf1(i,j,1)
            sy1=blk%surf2(i,j-1,1)
            sy2=blk%surf2(i,j,1)
            vol=blk%vol(i,j,1)
            u=blk%u0(i,j,1,:)
            xflux1=blk%xflux(i-1,j,1,1:5)
            xflux2=blk%xflux(i,j,1,1:5)
            yflux1=blk%yflux(i,j-1,1,1:5)
            yflux2=blk%yflux(i,j,1,1:5)
            blk%u(i,j,1,1:5)=u+(xflux1*sx1-xflux2*sx2+yflux1*sy1-yflux2*sy2)*dt/vol
        end do
    end do
    call convert_u_to_w_block(blk)
end subroutine vanleer_corrector

end module muscl
