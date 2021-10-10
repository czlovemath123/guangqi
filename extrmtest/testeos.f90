program testeos
use eosriemannsolver
use mathlib
implicit none
real(8), dimension(4) :: wl,wr,w
procedure(fun), pointer :: ptr
real(8) :: pstar,ustar,ws,rhostar(2),lambda(2),dxdt,  &
rhoscale,vscale,tempscale,vmax,rhol,rhor,vl,vr,pl,pr,gamml,gammr,  &
rho_temp,p_temp,t_temp,gamm,pressure,rho,cs,log10_rho,log10_temp,log10_p,  &
log10_cs,eg,converge,ws2min,ws2max,temperature,templ,tempr,log10_templ,  &
log10_tempr,log10_rhol,log10_rhor,log10_pl,log10_pr,mu,dt,t1,t2,  &
species(7),h_ratio,he_ratio
real(8), dimension(:), allocatable :: x,root,rholist,tlist,plist,temp1,temp2,temp3
integer(4) :: i,j,k,nrho,n,o,boundl(4),boundr(4),nv,n_sample
character(len=30) :: fmt1,fmt2,fmt3,subdir1,subdir2
character(len=4) :: mark1,mark2,mark3
character(len=5) :: str1,str2,str3,str4
real(8) :: test(4,6),test1(6),test2(6),test3(6),test4(6),dv(4),dv_const,dlogp
real(8), allocatable :: ustarl(:),ustarr(:),pstar_guess(:)
namelist /example/ test1,test2,test3,test4,h_ratio,he_ratio
namelist /exact/ test1,vl,vr,nv
#if     puphase==1
    !plot p-u curves
    open(unit=11,file="modules/eos_testsuit/problem.data",status="old",action="read")
    read(unit=11,nml=example)
    close(11)
    open(unit=22,file='modules/eos_testsuit/eos_exact/puphase.dat',status='replace',action='write')
    nrho=801
    n=2501
    o=nrho+n
    environment%table_dim=(/nrho,n/)
    environment%meta_table_dim=(/801,1401,1401/)
    write(str1,'(I4)')nrho
    write(str2,'(I4)')n
    write(str3,'(F4.2)')h_ratio
    write(str4,'(F4.2)')he_ratio
    str1=adjustl(str1)
    str2=adjustl(str2)
    str3=adjustl(str3)
    str4=adjustl(str4)
    subdir1='x'//trim(str3)//'_y'//trim(str4)//'/'
    subdir2='mn'//trim(str1)//'_'//trim(str2)
    call generate_eos_tables(subdir1,subdir2)
    rhol=test1(1)
    rhor=test1(4)
    templ=test1(3)
    tempr=test1(6)
    log10_templ=log10(templ)
    log10_tempr=log10(tempr)
    log10_rhol=log10(rhol)
    log10_rhor=log10(rhor)
    call interpolation_linear(log10_templ,log10_rhol,log10_pl,p_trho%xlist,p_trho%ylist,p_trho%t2d)
    call interpolation_linear(log10_tempr,log10_rhor,log10_pr,p_trho%xlist,p_trho%ylist,p_trho%t2d)
    pl=10**log10_pl
    pr=10**log10_pr
    gamml=1.667
    gammr=1.667
    wl=(/rhol,test1(2),pl,gamml/)
    wr=(/rhor,test1(5),pr,gammr/)
    call eos_riemann(wl,wr,pstar,ustar,rhostar,lambda)
    print *,pstar,ustar
    call eos_riemann_secant(wl,wr,pstar,ustar)
    print *,pstar,ustar
    n_sample=101
    dlogp=log10(max(wl(3),wr(3))/min(wl(3),wr(3)))/(n_sample-1)
    print *,'puphase: ','pmax=',max(wl(3),wr(3)),'pmin=',min(wl(3),wr(3)),dlogp
    allocate(ustarl(n_sample),ustarr(n_sample),pstar_guess(n_sample))
    do i=1,n_sample
        pstar_guess(i)=min(wl(3),wr(3))*10**((i-1)*dlogp)
        ustarl(i)=eos_exactrm_ustarl(pstar_guess(i),wl,wr)
        ustarr(i)=eos_exactrm_ustarr(pstar_guess(i),wl,wr)
        write(*,fmt='(5ES15.4E2)') pstar_guess(i),ustarl(i),ustarr(i)
        write(22,fmt='(5ES15.4E2)') pstar_guess(i),ustarl(i),ustarr(i)
    end do
    deallocate(ustarl,ustarr,pstar_guess)
    close(22)
#else
    !solve riemann problem
#if     ieos==1
    !use constant gamma gas eos table
    subdir1='constant_gamma/'
    environment%table_dim=(/400,601/)
    environment%meta_table_dim=(/501,701,701/)
    !constant gamma EoS table
#if     testid==1
    subdir2='gamma1.05'
    open(unit=23,file='modules/eos_testsuit/eos_exact/gamma1.dat',status='replace',action='write')
#elif   testid==2
    subdir2='gamma1.667'
    open(unit=23,file='modules/eos_testsuit/eos_exact/gamma2.dat',status='replace',action='write')
#endif
    open(unit=11,file='modules/extriemanntest/problem1.data',status='old',action='read')
    read(unit=11,nml=exact)
    close(11)
    dv_const=(vr-vl)/(nv-1)
    call generate_eos_tables(subdir1,subdir2)
    rhol=test1(1)
    rhor=test1(4)
    templ=test1(3)
    tempr=test1(6)
    log10_templ=log10(templ)
    log10_tempr=log10(tempr)
    log10_rhol=log10(rhol)
    log10_rhor=log10(rhor)
    call interpolation_linear(log10_templ,log10_rhol,log10_pl,p_trho%xlist,p_trho%ylist,p_trho%t2d)
    call interpolation_linear(log10_tempr,log10_rhor,log10_pr,p_trho%xlist,p_trho%ylist,p_trho%t2d)
    pl=10**log10_pl
    pr=10**log10_pr
    gamml=1.4
    gammr=1.4
    wl=(/rhol,test1(2),pl,gamml/)
    wr=(/rhor,test1(5),pr,gammr/)
    print *,wl
    print *,wr
    call eos_riemann(wl,wr,pstar,ustar,rhostar,lambda)
    print *,pstar
    do j=0,nv-1
        dxdt=vl+j*dv_const
        call eos_xstate(wl,wr,w,dxdt,pstar,ustar,rhostar,lambda)
        log10_p=log10(w(3))
        log10_rho=log10(w(1))
        call interpolation_linear(log10_p,log10_rho,log10_temp,t_prho%xlist,t_prho%ylist,t_prho%t2d)
        call interpolation_linear(log10_temp,log10_rho,gamm,gamma_trho%xlist,gamma_trho%ylist,gamma_trho%t2d)
        call interpolation_linear(log10_temp,log10_rho,mu,maw_trho%xlist,maw_trho%ylist,maw_trho%t2d)
        call interpolation_linear(log10_temp,log10_rho,eg,eg_trho%xlist,eg_trho%ylist,eg_trho%t2d)
        !this eg is the log10 of physical eg in cgs unit
        write(unit=23,fmt='(15ES18.6E2)') dxdt,w(1),w(2),w(3),10**log10_temp
        print *,dxdt/1e5,w(1),w(2)/1e5,w(3),10**log10_temp
    end do
    close(23)
#elif   ieos==2
    !use realistic gas eos table
    open(unit=11,file="modules/eos_testsuit/problem.data",status="old",action="read")
    open(unit=23,file='modules/eos_testsuit/eos_exact/eosreal1.dat',status='replace',action='write')
    open(unit=24,file='modules/eos_testsuit/eos_exact/eosreal2.dat',status='replace',action='write')
    open(unit=25,file='modules/eos_testsuit/eos_exact/eosreal3.dat',status='replace',action='write')
    open(unit=26,file='modules/eos_testsuit/eos_exact/eosreal4.dat',status='replace',action='write')
    read(unit=11,nml=example)
    close(11)
    test(1,:)=test1
    test(2,:)=test2
    test(3,:)=test3
    test(4,:)=test4
    nrho=401
    n=451
    o=nrho+n
    environment%table_dim=(/nrho,n/)
    environment%meta_table_dim=(/501,701,701/)
    write(str1,'(I4)')nrho
    write(str2,'(I4)')n
    write(str3,'(F4.2)')h_ratio
    write(str4,'(F4.2)')he_ratio
    str1=adjustl(str1)
    str2=adjustl(str2)
    str3=adjustl(str3)
    str4=adjustl(str4)
    subdir1='x'//trim(str3)//'_y'//trim(str4)//'/'
    subdir2='mn'//trim(str1)//'_'//trim(str2)
    call generate_eos_tables(subdir1,subdir2)
    boundl=(/-300,-200,-275,-320/)
    boundr=(/400,200,275,320/)
    dv=(/1e4,5e3,2e4,2.5e3/)
    do i=1,2
        rhol=test(i,1)
        rhor=test(i,4)
        templ=test(i,3)
        tempr=test(i,6)
        log10_templ=log10(templ)
        log10_tempr=log10(tempr)
        log10_rhol=log10(rhol)
        log10_rhor=log10(rhor)
        call interpolation_linear(log10_templ,log10_rhol,log10_pl,p_trho%xlist,p_trho%ylist,p_trho%t2d)
        call interpolation_linear(log10_tempr,log10_rhor,log10_pr,p_trho%xlist,p_trho%ylist,p_trho%t2d)
        pl=10**log10_pl
        pr=10**log10_pr
        gamml=1.4
        gammr=1.4
        wl=(/rhol,test(i,2),pl,gamml/)
        wr=(/rhor,test(i,5),pr,gammr/)
        print *,wl
        print *,wr
        call eos_riemann(wl,wr,pstar,ustar,rhostar,lambda)
        call cpu_time(t1)
        k=22+i
        write(unit=k,fmt='(19A16)') 'dxdt','rho','vx','p','temperature','mu','eg','Gamma','2nH2/nHtot','nHI/nHtot',  &
            'nHII/nHtot','nHeI/nHetot','nHeII/nHetot','nHeIII/nHetot','ne/nHtot'
        do j=boundl(i),boundr(i)
            dxdt=j*dv(i)
            call eos_xstate(wl,wr,w,dxdt,pstar,ustar,rhostar,lambda)
            log10_p=log10(w(3))
            log10_rho=log10(w(1))
            call interpolation_linear(log10_p,log10_rho,log10_temp,t_prho%xlist,t_prho%ylist,t_prho%t2d)
            call interpolation_linear(log10_temp,log10_rho,gamm,gamma_trho%xlist,gamma_trho%ylist,gamma_trho%t2d)
            call interpolation_linear(log10_temp,log10_rho,mu,maw_trho%xlist,maw_trho%ylist,maw_trho%t2d)
            call interpolation_linear(log10_temp,log10_rho,eg,eg_trho%xlist,eg_trho%ylist,eg_trho%t2d)
            call cal_species_h_he_indivi(h_ratio,he_ratio,mu,species)
            !this eg is the log10 of physical eg in cgs unit
            write(unit=k,fmt='(15ES16.6E2)') dxdt,w(1),w(2),w(3),10**log10_temp,mu,eg,gamm,species
            print *,dxdt/1e5,w(1),w(2)/1e5,w(3),10**log10_temp,mu,eg
        end do
        call cpu_time(t2)
        dt=t2-t1
        print *,'eos exact riemann solver, test',i,'t=',dt
        close(k)
    end do
#endif
#endif
end program testeos
