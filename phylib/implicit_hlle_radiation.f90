module implicit_hlle_radiation
use phylib
use mathlib
use datastructure
use eos
use radiation_common_functions
implicit none

contains

!subroutine radiation_hlle_implicit(info,dt)
!    !solve Erad and Fradx implicitly
!    type(infodef) :: info
!    real(8) :: theta(6),phi(6),b(nx2),a(nx2,nx2),ab(ktotal,nx2),alpha,dt,dx,r,chi,fl,fm,fr,  &
!        rthetaphi(3),gradvf,s,temp,r_ratio,qrad,Erad,egv,egv_predict,temp_predict,rho,  &
!        Fradxp,Fradxn,dxyz_block(3),pos_block(3)
!    real(8), dimension(:,:,:,:), allocatable :: u
!    integer :: i,j,k,ijk(3),err,nrhs,ldab,ldb,ipiv(nx2)
!    dx=dxyz(1)
!    alpha=dt/dx
!    b=0d0
!    if (nd==1) then
!        call allocate_cell_data(u,3)
!        u=info%w(xlb:xub,1:1,1:1,1:3)
!    else if (nd==2) then
!        call allocate_cell_data(u,3)
!        u=info%w(xlb:xub,ylb:yub,1:1,1:3)
!    end if
!    call gradient(u,info%gradv,dxyz(1))
!    do i=1,nx
!        if (info%Erad(i,1,1)<0) then
!            !info%Erad(i,1,1)=1d-3
!            write(*,'(10ES18.5E3)') info%Erad
!            write(*,'(10ES18.5E3)') info%Fradx
!            print *,'density'
!            write(*,'(10ES18.5E3)') info%w(:,1,1,1)
!            print *,'temperature'
!            write(*,'(10ES18.5E3)') info%temp
!            print *,'velocity'
!            write(*,'(10ES18.5E3)') info%w(:,1,1,2)
!            print *,'opacity'
!            write(*,'(10ES18.5E3)') info%opacity
!            print *,'hlle implicit1',alpha
!            stop
!        end if
!    end do
!    if (nd==1) then
!        info%Erad_new=info%Erad
!        info%Fradx_new=info%Fradx
!        do i=1,nx
!            chi=info%chi(i,1,1)
!            temp=info%temp(i,1,1)
!            Erad=info%Erad(i,1,1)
!            Fradxn=info%Fradx(i-1,1,1)
!            Fradxp=info%Fradx(i+1,1,1)
!            qrad=chi*c_light*(a_rad*temp**4-Erad)!-c_light*(Fradxp-Fradxn)/dx/2
!            egv=info%egv(i,1,1)
!            egv_predict=egv-qrad*dt
!            rho=info%w(i,1,1,1)
!            temp_predict=calculate_temp_from_egv_rho(egv_predict,rho)
!            !print *,temp,temp_predict,qrad,a_rad*temp**4/Erad
!            !s=4d0*pi*(planck_function(temp)+planck_function(temp_predict))/2*chi*dt
!            s=4d0*pi*planck_function(temp)*chi*dt
!            if (i==1) then
!                if (rad_bound_type(1)==transmissive) then
!                    b(1)=info%Erad(1,1,1)+s
!                    b(2)=info%Fradx(1,1,1)
!                else if (rad_bound_type(1)==reflective) then
!                    !b(1)=info%Erad(1,1,1)+s
!                    !b(2)=info%Fradx(1,1,1)
!                else if (rad_bound_type(1)==specified) then
!                    !theta(1)=theta1(alpha,info%feddington(0,1,1),info%feddington(1,1,1))
!                    !theta(4)=theta4(alpha,info%feddington(0,1,1),info%feddington(1,1,1))
!                    !phi(1)=phi1(alpha,info%feddington(0,1,1),info%feddington(1,1,1))
!                    !phi(4)=phi4(alpha,info%feddington(0,1,1),info%feddington(1,1,1))
!                    !b(1)=info%Erad(1,1,1)-theta(1)*info%Erad(0,1,1)-theta(4)*info%Fradx(0,1,1)+s
!                    !b(2)=info%Fradx(0,1,1)-phi(1)*info%Erad(0,1,1)-phi(4)*info%Fradx(0,1,1)
!                end if
!            else if (i==nx) then
!                if (rad_bound_type(2)==transmissive) then
!                    b(nx2-1)=info%Erad(nx,1,1)+s
!                    b(nx2)=info%Fradx(nx,1,1)
!                else if (rad_bound_type(2)==reflective) then
!                    !b(nx2-1)=info%Erad(nx,1,1)+s
!                    !b(nx2)=info%Fradx(nx,1,1)
!                else if (rad_bound_type(2)==specified) then
!                    !theta(3)=theta3(alpha,info%feddington(nx,1,1),info%feddington(nx+1,1,1))
!                    !theta(6)=theta6(alpha,info%feddington(nx,1,1),info%feddington(nx+1,1,1))
!                    !phi(3)=phi3(alpha,info%feddington(nx,1,1),info%feddington(nx+1,1,1))
!                    !phi(6)=phi6(alpha,info%feddington(nx,1,1),info%feddington(nx+1,1,1))
!                    !b(nx2-1)=info%Erad(nx,1,1)-theta(3)*info%Erad(nx+1,1,1)-theta(6)*info%Fradx(nx+1,1,1)+s
!                    !b(nx2)=info%Fradx(nx,1,1)-phi(3)*info%Erad(nx+1,1,1)-phi(6)*info%Fradx(nx+1,1,1)
!                end if
!            else
!                b(2*(i-1)+1)=info%Erad(i,1,1)+s
!                b(2*(i-1)+2)=info%Fradx(i,1,1)
!            end if
!        end do
!        a=0d0
!        do i=1,nx
!            ijk=(/i,1,1/)
!            call ijk_to_coords(ijk,dxyz_block,pos_block,rthetaphi)
!            chi=info%chi(i,1,1)
!            if (i==1) then
!                fl=info%feddington(0,1,1)
!                fm=info%feddington(1,1,1)
!                fr=info%feddington(2,1,1)
!                if (rad_bound_type(1)==transmissive) then
!                    theta(1)=theta1(alpha,fl,fm)
!                    theta(2)=theta2(alpha,fl,fm)
!                    theta(3)=theta3(alpha,fl,fm,fr,chi,dt)
!                    theta(4)=theta4(alpha,fl,fm,fr)
!                    theta(5)=theta5(alpha,fm,fr)
!                    theta(6)=theta6(alpha,fm,fr)
!                    phi(1)=phi1(alpha,fl,fm)
!                    phi(2)=phi2(alpha,fl,fm)
!                    phi(3)=phi3(alpha,fl,fm,fr)
!                    phi(4)=phi4(alpha,fl,fm,fr,chi,dt)
!                    phi(5)=phi5(alpha,fm,fr)
!                    phi(6)=phi6(alpha,fm,fr)
!                    a(1,1)=theta(1)+theta(3)
!                    a(1,2)=theta(2)+theta(4)
!                    a(1,3)=theta(5)
!                    a(1,4)=theta(6)
!                    a(2,1)=phi(1)+phi(3)
!                    a(2,2)=phi(2)+phi(4)
!                    a(2,3)=phi(5)
!                    a(2,4)=phi(6)
!                else if (rad_bound_type(1)==reflective) then
!                    !theta(1)=theta1(alpha,fl,fm)
!                    !theta(2)=theta2(alpha,fl,fm,fr,gradvf,chi,dt)
!                    !theta(3)=theta3(alpha,fm,fr)
!                    !theta(4)=theta4(alpha,fl,fm)
!                    !theta(5)=theta5(alpha,fl,fm,fr,dt,r)
!                    !theta(6)=theta6(alpha,fm,fr)
!                    !phi(1)=phi1(alpha,fl,fm)
!                    !phi(2)=phi2(alpha,fl,fm,fr,dt,r)
!                    !phi(3)=phi3(alpha,fm,fr)
!                    !phi(4)=phi4(alpha,fl,fm)
!                    !phi(5)=phi5(alpha,fl,fm,fr,chi,dt)
!                    !phi(6)=phi6(alpha,fm,fr)
!                    !a(1,1)=theta(2)+theta(1)
!                    !a(1,2)=theta(5)-theta(4)
!                    !a(1,3)=theta(3)
!                    !a(1,4)=theta(6)
!                    !a(2,1)=phi(2)+phi(1)
!                    !a(2,2)=phi(5)-phi(4)
!                    !a(2,3)=phi(3)
!                    !a(2,4)=phi(6)
!                else if (rad_bound_type(1)==specified) then
!                    !theta(2)=theta2(alpha,fl,fm,fr,gradvf,chi,dt)
!                    !theta(3)=theta3(alpha,fm,fr)
!                    !theta(5)=theta5(alpha,fl,fm,fr,dt,r)
!                    !theta(6)=theta6(alpha,fm,fr)
!                    !phi(2)=phi2(alpha,fl,fm,fr,dt,r)
!                    !phi(3)=phi3(alpha,fm,fr)
!                    !phi(5)=phi5(alpha,fl,fm,fr,chi,dt)
!                    !phi(6)=phi6(alpha,fm,fr)
!                    !a(1,1)=theta(2)
!                    !a(1,2)=theta(5)
!                    !a(1,3)=theta(3)
!                    !a(1,4)=theta(6)
!                    !a(2,1)=phi(2)
!                    !a(2,2)=phi(5)
!                    !a(2,3)=phi(3)
!                    !a(2,4)=phi(6)
!                end if
!            else if (i==nx) then
!                fl=info%feddington(nx-1,1,1)
!                fm=info%feddington(nx,1,1)
!                fr=info%feddington(nx+1,1,1)
!                if (rad_bound_type(2)==transmissive) then
!                    theta(1)=theta1(alpha,fl,fm)
!                    theta(2)=theta2(alpha,fl,fm)
!                    theta(3)=theta3(alpha,fl,fm,fr,chi,dt)
!                    theta(4)=theta4(alpha,fl,fm,fr)
!                    theta(5)=theta5(alpha,fm,fr)
!                    theta(6)=theta6(alpha,fm,fr)
!                    phi(1)=phi1(alpha,fl,fm)
!                    phi(2)=phi2(alpha,fl,fm)
!                    phi(3)=phi3(alpha,fl,fm,fr)
!                    phi(4)=phi4(alpha,fl,fm,fr,chi,dt)
!                    phi(5)=phi5(alpha,fm,fr)
!                    phi(6)=phi6(alpha,fm,fr)
!                    a(nx2-1,nx2-3)=theta(1)
!                    a(nx2-1,nx2-2)=theta(2)
!                    a(nx2-1,nx2-1)=theta(3)+theta(5)
!                    a(nx2-1,nx2)=theta(4)+theta(6)
!                    a(nx2,nx2-3)=phi(1)
!                    a(nx2,nx2-2)=phi(2)
!                    a(nx2,nx2-1)=phi(3)+phi(5)
!                    a(nx2,nx2)=phi(4)+phi(6)
!                else if (rad_bound_type(2)==reflective) then
!                else if (rad_bound_type(2)==extrapolate) then
!                    !theta(1)=theta1(alpha,fl,fm)
!                    !theta(2)=theta2(alpha,fl,fm,fr,gradvf,chi,dt)
!                    !theta(3)=theta3(alpha,fm,fr)
!                    !theta(4)=theta4(alpha,fl,fm)
!                    !theta(5)=theta5(alpha,fl,fm,fr,dt,r)
!                    !theta(6)=theta6(alpha,fm,fr)
!                    !phi(1)=phi1(alpha,fl,fm)
!                    !phi(2)=phi2(alpha,fl,fm,fr,dt,r)
!                    !phi(3)=phi3(alpha,fm,fr)
!                    !phi(4)=phi4(alpha,fl,fm)
!                    !phi(5)=phi5(alpha,fl,fm,fr,chi,dt)
!                    !phi(6)=phi6(alpha,fm,fr)
!                    !a(nx2-1,nx2-3)=theta(1)-theta(3)
!                    !a(nx2-1,nx2-2)=theta(4)-theta(6)
!                    !a(nx2-1,nx2-1)=theta(2)+2*theta(3)
!                    !a(nx2-1,nx2)=theta(5)+2*theta(6)
!                    !a(nx2,nx2-3)=phi(1)-phi(3)
!                    !a(nx2,nx2-2)=phi(4)-phi(6)
!                    !a(nx2,nx2-1)=phi(2)+2*phi(3)
!                    !a(nx2,nx2)=phi(5)+2*phi(6)
!                else if (rad_bound_type(2)==specified) then
!                    !theta(1)=theta1(alpha,fl,fm)
!                    !theta(2)=theta2(alpha,fl,fm,fr,gradvf,chi,dt)
!                    !theta(4)=theta4(alpha,fl,fm)
!                    !theta(5)=theta5(alpha,fl,fm,fr,dt,r)
!                    !phi(1)=phi1(alpha,fl,fm)
!                    !phi(2)=phi2(alpha,fl,fm,fr,dt,r)
!                    !phi(4)=phi4(alpha,fl,fm)
!                    !phi(5)=phi5(alpha,fl,fm,fr,chi,dt)
!                    !a(nx2-1,nx2-3)=theta(1)
!                    !a(nx2-1,nx2-2)=theta(4)
!                    !a(nx2-1,nx2-1)=theta(2)
!                    !a(nx2-1,nx2)=theta(5)
!                    !a(nx2,nx2-3)=phi(1)
!                    !a(nx2,nx2-2)=phi(4)
!                    !a(nx2,nx2-1)=phi(2)
!                    !a(nx2,nx2)=phi(5)
!                end if
!            else
!                fl=info%feddington(i-1,1,1)
!                fm=info%feddington(i,1,1)
!                fr=info%feddington(i+1,1,1)
!                theta(1)=theta1(alpha,fl,fm)
!                theta(2)=theta2(alpha,fl,fm)
!                theta(3)=theta3(alpha,fl,fm,fr,chi,dt)
!                theta(4)=theta4(alpha,fl,fm,fr)
!                theta(5)=theta5(alpha,fm,fr)
!                theta(6)=theta6(alpha,fm,fr)
!                phi(1)=phi1(alpha,fl,fm)
!                phi(2)=phi2(alpha,fl,fm)
!                phi(3)=phi3(alpha,fl,fm,fr)
!                phi(4)=phi4(alpha,fl,fm,fr,chi,dt)
!                phi(5)=phi5(alpha,fm,fr)
!                phi(6)=phi6(alpha,fm,fr)
!                k=2*(i-1)+1
!                a(k,k-2)=theta(1)
!                a(k,k-1)=theta(2)
!                a(k,k)=theta(3)
!                a(k,k+1)=theta(4)
!                a(k,k+2)=theta(5)
!                a(k,k+3)=theta(6)
!                k=2*(i-1)+2
!                a(k,k-3)=phi(1)
!                a(k,k-2)=phi(2)
!                a(k,k-1)=phi(3)
!                a(k,k)=phi(4)
!                a(k,k+1)=phi(5)
!                a(k,k+2)=phi(6)
!            end if
!        end do
!        ab=0d0
!        do i=1,nx2
!            do j=1,nx2
!                if (i>=max(1,j-ku_hlle).and.i<=min(nx2,j+kl_hlle)) then
!                    ab(kl_hlle+ku_hlle+1+i-j,j)=a(i,j)
!                end if
!            end do
!        end do
!        ldb=nx2
!        nrhs=1
!        ldab=ktotal
!        !print *,'a matrix'
!        !write(*,'(20ES10.2E2)') a
!        !print *,'band storage'
!        !write(*,'(10ES10.2E2)') ab
!        !print *,'rhs'
!        !write(*,'(10ES12.3E2)') b
!                !write(*,'(10ES18.5E3)') info%Erad
!        call dgbsv(nx2,kl_hlle,ku_hlle,nrhs,ab,ldab,ipiv,b,ldb,err)
!        if (err/=0) then
!            print *,err,'radiation hlle implicit wrong'
!            stop
!        end if
!        do i=1,nx
!            if (b(2*(i-1)+1)<0) then
!                !info%Erad(i,1,1)=1d-5
!                !print *,i
!                write(*,'(10ES18.5E3)') info%Erad
!                write(*,'(10ES18.5E3)') info%Fradx
!                print *,'density'
!                write(*,'(10ES18.5E3)') info%w(:,1,1,1)
!                print *,'temperature'
!                write(*,'(10ES18.5E3)') info%temp
!                print *,'velocity'
!                write(*,'(10ES18.5E3)') info%w(:,1,1,2)
!                print *,'opacity'
!                write(*,'(10ES18.5E3)') info%opacity
!                print *,'hlle implicit2',alpha
!                stop
!            end if
!        end do
!        do i=1,nx
!            !info%Erad_new(i,1,1)=b(2*(i-1)+1)
!            !info%Fradx_new(i,1,1)=b(2*(i-1)+2)
!            info%Erad(i,1,1)=b(2*(i-1)+1)
!            info%Fradx(i,1,1)=b(2*(i-1)+2)
!        end do
!        !print *,info%Erad
!    else
!        !2d, not implemented yet
!    end if
!    deallocate(u)
!end subroutine radiation_hlle_implicit
!
!subroutine hlle_radiation_momentum_transfer(info,dt)
!    !radiation flux is calculated at the interfaces
!    !flux value at cell center is interpolated
!    type(infodef) :: info
!    real(8) :: dt,dx,chi,chi_x,chi_y,lambda_x,lambda_y,dx_coefficient,dy_coefficient,du(5),rho,dek,  &
!        v_old(3),v_new(3)
!    integer :: i,j
!    do i=1,nx
!        if (info%temp(i,1,1)<0) then
!            print *,i,'momentum transfer1'
!            print *,info%temp
!            stop
!        end if
!    end do
!    dx=dxyz(1)
!    if (nd==1) then
!        do i=xlb,xub
!            rho=info%w(i,1,1,irho)
!            v_old=info%w(i,1,1,ivx:ivz)
!            du=zero
!            chi=info%chi(i,1,1)
!            du(2)=(info%Fradx(i,1,1)+info%Fradx_new(i,1,1))/2*chi*dt
!            v_new=v_old+du(2:4)/rho
!            if (abs(v_new(1))>3d7) then
!                !v_new(1)=1d7
!                !du(2)=rho*(v_new(1)-v_old(1))
!                print *,'too fast or momentum stiff'
!                write(*,'(10ES18.5E3)') info%opacity
!                print *,i,info%opacity(i,1,1),v_old(1),v_new(1),info%Fradx(i,1,1)
!                stop
!            end if
!            dek=half*rho*(v_new(1)**2-v_old(1)**2)
!            du(5)=dek
!            info%u(i,1,1,1:5)=info%u(i,1,1,1:5)+du
!        end do
!        info%Erad=info%Erad_new
!        info%Fradx=info%Fradx_new
!        info%Fradx_total=info%Fradx_total+info%Fradx_new*dt
!    else if (nd==2) then
!        do j=1,ny
!            do i=1,nx
!                rho=info%w(i,j,1,irho)
!                v_old=info%w(i,j,1,ivx:ivz)
!                du=zero
!                chi=info%chi(i,j,1)
!                du(2)=info%Fradx(i,j,1)*chi*dt
!                du(3)=info%Frady(i,j,1)*chi*dt
!                v_new=v_old+du(2:4)/rho
!                dek=half*rho*(dot_product(v_new,v_new)-dot_product(v_old,v_old))
!                du(5)=dek
!                info%u(i,j,1,1:5)=info%u(i,j,1,1:5)+du
!            end do
!        end do
!    else
!    end if
!    call convert_u_to_w_block(info%u,info%w,info%temp,info%egv)
!    do i=1,nx
!        if (info%temp(i,1,1)<0) then
!            print *,i,'momentum transfer2'
!            print *,info%temp
!            print *,info%w(:,1,1,ivx)
!            stop
!        end if
!    end do
!end subroutine hlle_radiation_momentum_transfer
!
!subroutine hlle_radiation_matter_interaction(info,dt)
!    type(infodef) :: info
!    real(8) :: chi,dx,egn,eradn,tempn,dt,egnext,rho,eradnext,ek,Fradxp,Fradxn,temp_predict,egv_predict,qrad
!    integer :: i,j
!    dx=dxyz(1)
!    do i=1,nx
!        if (info%temp(i,1,1)<0) then
!            print *,i,'radiation matter1'
!            print *,info%temp
!            stop
!        end if
!    end do
!    if (nd==1) then
!        do i=1,nx
!            chi=info%chi(i,1,1)
!            if (chi/=0) then
!                egn=info%egv(i,1,1)
!                rho=info%w(i,1,1,irho)
!                eradn=(info%Erad(i,1,1)+info%Erad(i,1,1))/2
!                !egnext=radiation_matter_egnext(chi,egn,eradn,rho,dt)
!                tempn=info%temp(i,1,1)
!                !Fradxn=info%Fradx(i-1,1,1)
!                !Fradxp=info%Fradx(i+1,1,1)
!                !qrad=chi*c_light*(a_rad*tempn**4-eradn)!-c_light*(Fradxp-Fradxn)/dx/2
!                !egv_predict=egn-qrad*dt
!                !temp_predict=calculate_temp_from_egv_rho(egv_predict,rho)
!                egnext=egn-chi*c_light*(a_rad*tempn**4-eradn)*dt!-c_light*(Fradxp-Fradxn)/2/dx
!                ek=info%u(i,1,1,5)-egn
!                !print *,egnext
!                info%u(i,1,1,5)=ek+egnext
!            end if
!        end do
!    else if (nd==2) then
!        do j=1,ny
!            do i=1,nx
!                chi=info%chi(i,j,1)
!                egn=info%egv(i,j,1)
!                rho=info%w(i,j,1,irho)
!                eradn=info%Erad(i,j,1)
!                egnext=hlle_radiation_matter_egnext(chi,egn,eradn,rho,dt)
!                ek=info%u(i,j,1,5)-egn
!                info%u(i,j,1,5)=ek+egnext
!            end do
!        end do
!    else
!    end if
!    call convert_u_to_w_block(info%u,info%w,info%temp,info%egv)
!    !print *,'2',info%temp
!    !stop
!    do i=1,nx
!        if (info%temp(i,1,1)<0.or.info%w(i,1,1,irho)<0.or.info%w(i,1,1,5)<0.or.info%egv(i,1,1)<0.or.isnan(info%temp(i,1,1))) then
!            print *,i,'radiation matter2'
!            print *,info%temp
!            stop
!        end if
!        if (info%Erad(i,1,1)<0) then
!            print *,'radiation matter3'
!        end if
!    end do
!end subroutine hlle_radiation_matter_interaction
!
!subroutine radiation_ddt_estimate(info,dt_remain,ddt)
!    type(infodef) :: info
!    real(8) :: dt_remain,ddt,dx,egv,qrad,dt_cool,temp,chi,Erad,vx,dt_mom,Fradx,a,Fradxn,Fradxp
!    integer :: i,j,imin
!    ddt=dt_remain
!    dx=dxyz(1)
!    imin=0
!    if (nd==1) then
!        do i=1,nx
!            egv=info%egv(i,1,1)
!            chi=info%chi(i,1,1)
!            temp=info%temp(i,1,1)
!            Erad=info%Erad(i,1,1)
!            Fradx=info%Fradx(i,1,1)
!            vx=info%w(i,1,1,1)
!            Fradxn=info%Fradx(i-1,1,1)
!            Fradxp=info%Fradx(i+1,1,1)
!            qrad=chi*c_light*(a_rad*temp**4-Erad)!+c_light*(Fradxp-Fradxn)/2/dx
!            a=chi*Fradx
!            dt_cool=abs(egv/qrad/10)
!            dt_mom=max(abs(vx),1d6)/max(abs(a),1d0)
!            !print *,i,dt_cool,dt_mom
!            if (dt_cool<ddt.or.dt_mom<ddt) then
!                ddt=dt_cool
!                imin=i
!            end if
!        end do
!        !print *,imin,(info%Erad(imin,1,1)/a_rad)**0.25d0,info%temp(imin,1,1),dt_remain
!    else
!    end if
!end subroutine radiation_ddt_estimate
!
!function hlle_radiation_matter_egnext(chi,egn,erad,rho,dt)
!        !use implicit method to calculate the gas internal energy of the next timestep
!    real(8) :: egn,erad,chi,dt,rho
!    real(8) :: hlle_radiation_matter_egnext
!    real(8) :: a(5),a1,a2
!    complex(8) :: z(4)
!    a1=chi*a_rad*c_light*dt*((gamma_gas-one)*maw*amu/rho/kb)**4
!    a2=c_light*chi*erad*dt+egn
!    a(1)=-a2/a1
!    a(2)=one/a1
!    a(3)=zero
!    a(4)=zero
!    a(5)=one
!    call quarticroots(a,z)
!    hlle_radiation_matter_egnext=dble(z(1))
!end function hlle_radiation_matter_egnext
!
!subroutine hlle_calculate_feddington(info)
!    type(infodef) :: info
!    real(8) :: El,Er,Erad,grad_Erad,dx,chi,r,lambda
!    integer :: i
!    dx=dxyz(1)
!    do i=1,nx
!        info%feddington=1d0
!    end do
!end subroutine hlle_calculate_feddington
!
!function theta1(alpha,fl,fr)
!    real(8) :: alpha,fl,fr,theta1
!    theta1=-alpha*c_light*sqrt(fl*fr)/(sqrt(fl)+sqrt(fr))
!end function theta1
!
!function theta2(alpha,fl,fr)
!    real(8) :: alpha,fl,fr,theta2
!    theta2=-alpha*c_light*sqrt(fr)/(sqrt(fl)+sqrt(fr))
!end function theta2
!
!function theta3(alpha,fl,fm,fr,chi,dt)
!    real(8) :: alpha,fl,fm,fr,theta3,chi,dt
!    theta3=1d0+alpha*c_light*(sqrt(fl*fm)/(sqrt(fl)+sqrt(fm))+sqrt(fm*fr)/(sqrt(fm)+sqrt(fr)))  &
!        +chi*c_light*dt
!end function theta3
!
!function theta4(alpha,fl,fm,fr)
!    real(8) :: alpha,fl,fm,fr,theta4
!    theta4=alpha*c_light*(sqrt(fr)/(sqrt(fm)+sqrt(fr))-sqrt(fl)/(sqrt(fl)+sqrt(fm)))
!end function theta4
!
!function theta5(alpha,fl,fr)
!    real(8) :: alpha,fl,fr,theta5
!    theta5=-alpha*c_light*sqrt(fl*fr)/(sqrt(fl)+sqrt(fr))
!end function theta5
!
!function theta6(alpha,fl,fr)
!    real(8) :: alpha,fl,fr,theta6
!    theta6=alpha*c_light*sqrt(fl)/(sqrt(fl)+sqrt(fr))
!end function theta6
!
!function phi1(alpha,fl,fr)
!    real(8) :: alpha,fl,fr,phi1
!    phi1=-alpha*c_light*fl*sqrt(fr)/(sqrt(fl)+sqrt(fr))
!end function phi1
!
!function phi2(alpha,fl,fr)
!    real(8) :: alpha,fl,fr,phi2
!    phi2=-alpha*c_light*sqrt(fl*fr)/(sqrt(fl)+sqrt(fr))
!end function phi2
!
!function phi3(alpha,fl,fm,fr)
!    real(8) :: alpha,fl,fm,fr,phi3
!    phi3=alpha*c_light*(fm*sqrt(fr)/(sqrt(fm)+sqrt(fr))-fm*sqrt(fl)/(sqrt(fl)+sqrt(fm)))
!end function phi3
!
!function phi4(alpha,fl,fm,fr,chi,dt)
!    real(8) :: alpha,fl,fm,fr,phi4,chi,dt
!    phi4=1d0+alpha*c_light*(sqrt(fl*fm)/(sqrt(fl)+sqrt(fm))+sqrt(fm*fr)/(sqrt(fm)+sqrt(fr)))+chi*c_light*dt
!end function phi4
!
!function phi5(alpha,fl,fr)
!    real(8) :: alpha,fl,fr,phi5
!    phi5=alpha*c_light*fr*sqrt(fl)/(sqrt(fl)+sqrt(fr))
!end function phi5
!
!function phi6(alpha,fl,fr)
!    real(8) :: alpha,fl,fr,phi6
!    phi6=-alpha*c_light*sqrt(fl*fr)/(sqrt(fl)+sqrt(fr))
!end function phi6
!
!!not successful, the predictor-corrector scheme has huge problem
!
!subroutine stiff_source_picard_iteration(info,hydroflux)
!    !calculate info%u of the next timestep with Picard iteration
!    !radsource is in source form
!    type(infodef) :: info
!    real(8), dimension(:,:,:,:), allocatable :: hydroflux,u_next
!    real(8) :: temp,chi,rho,momx,e_total,c,s1,s2,s3,dt,Erad,egv,kappa
!    real(8), dimension(5) :: dudt_hydro,radsource,dudt_rad,dudt,u0,u_estimate,u_epsilon,u_corrected
!    integer :: i,j
!    dt=time_sys%dt_hydro
!    call allocate_cell_data(u_next,5)
!    if (nd==1) then
!        u_next=info%u
!        do i=xlb+1,xub-1
!            rho=info%u(i,1,1,1)
!            momx=info%u(i,1,1,2)
!            e_total=info%u(i,1,1,5)
!            temp=info%temp(i,1,1)
!            chi=info%chi(i,1,1)
!            kappa=info%opacity(i,1,1)
!            Erad=info%Erad(i,1,1)
!            radsource=zero
!            !radsource(5)=chi*c_light*(Erad-a_rad*temp**4)
!            radsource(5)=-c_light*(info%Fradx(i+1,1,1)-info%Fradx(i-1,1,1))/dxyz(1)/2
!            !print *,c_light*(info%Fradx(i+1,1,1)-info%Fradx(i-1,1,1))/dxyz(1)
!            c=4d0*chi*c_light*a_rad*temp**3*(gamma_gas-one)*maw*amu/kb
!            c=zero
!            s1=c*(-e_total/rho**2+momx**2/rho**3)
!            s2=-c*momx/rho**2
!            s3=c/rho
!            !print *,kappa*c_light*(Erad-a_rad*temp**4),s1,s2,s3
!            u0=info%u(i,1,1,1:5)
!            dudt_hydro=(hydroflux(i-1,1,1,1:5)-hydroflux(i,1,1,1:5))/dxyz(1)
!            dudt_rad=radsource
!            dudt=dudt_hydro+dudt_rad
!            !print *,dudt_hydro(5),dudt_rad(5),s3*dt
!            call stiff_source_picard_iteration_estimate(u0,s1,s2,s3,dudt,dt,u_estimate)
!            !call stiff_source_picard_iteration_error(Erad,chi,temp,u0,u_estimate,dudt_hydro,dt,u_epsilon)
!            !call stiff_source_picard_iteration_correction(u_estimate,u_epsilon,dt,chi,Erad,u_corrected)
!            u_next(i,1,1,1:5)=u_estimate!u_corrected
!        end do
!        info%u(xlb+1:xub-1,1,1,1:5)=u_next(xlb+1:xub-1,1,1,1:5)
!    else if (nd==2) then
!    else
!    end if
!    deallocate(u_next)
!end subroutine stiff_source_picard_iteration
!
!subroutine stiff_source_picard_iteration_estimate(u0,s1,s2,s3,dudt,dt,u_estimate)
!    !u0 the initial value
!    real(8) :: s1,s2,s3,dt
!    real(8), dimension(5) :: u_estimate,u0,dudt,du
!    if (nd==1) then
!        du=dudt*dt
!        u_estimate(1)=u0(1)+du(1)
!        u_estimate(2)=u0(2)+du(2)
!        u_estimate(3:4)=0d0
!        u_estimate(5)=u0(5)+(-s1*dt*du(1)-s2*dt*du(2)+du(5))/(1d0+s3*dt)
!        !print *,s3,s1/s3*du(1),s2/s3*du(2)
!        !print *,u_estimate(5)-half*u_estimate(2)**2/u_estimate(1),u0(5)-half*u0(2)**2/u0(1)
!        if ((u_estimate(5)-half*u_estimate(2)**2/u_estimate(1))<0) then
!            print *,'u_estimate negative'
!            print *,u0,(-s1*dt*du(1)-s2*dt*du(2)+du(5))/(1d0+s3*dt)
!            print *,s1*dt*du(1),s2*dt*du(2),du(5),1d0+s3*dt
!            !stop
!        end if
!        !print *,s1*dt*du(1),s2*dt*du(2),du(5),1d0-s3*dt
!    else
!    end if
!end subroutine stiff_source_picard_iteration_estimate
!
!subroutine stiff_source_picard_iteration_error(Erad,chi,temp1,u0,u_estimate,dudt_hydro,dt,u_epsilon)
!    real(8), dimension(5) :: u0,u_estimate,radsource1,radsource2,dudt_hydro,u_epsilon
!    real(8) :: Erad,chi,temp1,rho,egv,temp2,dt
!    radsource1=zero
!    radsource1(5)=chi*c_light*(Erad-a_rad*temp1**4)
!    rho=u_estimate(1)
!    egv=u_estimate(5)-half*u_estimate(2)**2/rho
!    temp2=calculate_temp_from_egv_rho(egv,rho)
!    !temp2=temp1
!    radsource1=zero
!    radsource2=zero
!    !radsource2(5)=chi*u_estimate(1)/u0(1)*c_light*(Erad-a_rad*temp2**4)
!    u_epsilon=u0-u_estimate+dt/2*(radsource1+radsource2)+dt*dudt_hydro
!end subroutine stiff_source_picard_iteration_error
!
!subroutine stiff_source_picard_iteration_correction(u_estimate,u_epsilon,dt,chi,Erad,u_corrected)
!    real(8), dimension(5) :: u_estimate,u_epsilon,u_corrected,delta
!    real(8) :: dt,chi,rho,momx,e_total,egv,c,temp,s1,s2,s3,Erad
!    rho=u_estimate(1)
!    momx=u_estimate(2)
!    e_total=u_estimate(5)
!    egv=e_total-half*momx**2/rho
!    temp=calculate_temp_from_egv_rho(egv,rho)
!    c=4d0*chi*c_light*a_rad*temp**3*(gamma_gas-one)*maw*amu/kb
!    c=zero
!    s1=c*(-e_total/rho**2+momx**2/rho**3)
!    s2=-c*momx/rho**2
!    s3=c/rho
!    delta=zero
!    delta(1)=u_epsilon(1)
!    delta(2)=u_epsilon(2)
!    delta(5)=(-s1*dt*u_epsilon(1)-s2*dt*u_epsilon(2)+u_epsilon(5))/(1d0+s3*dt)
!    u_corrected=u_estimate+delta
!end subroutine stiff_source_picard_iteration_correction

end module implicit_hlle_radiation
