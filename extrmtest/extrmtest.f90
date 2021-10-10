program extrmtest
use exactriemann
use phylib
real(8) :: test1(6)
real(8) :: xl,xr,yl,yr,zl,zr,dx,temperature,pl,pr,tl,tr,mu,lscale
integer(8) :: boundl,boundr
real(8), dimension(4) :: wl,wr,w
real(8) :: x,tolerance,root,ustar,t,ndensity,vl,vr,dlogp,gamma_law_eos
real(8), allocatable :: ustarl(:),ustarr(:),pstar_guess(:)
integer :: i,j,unitnumber,nv,n_sample
logical :: lriemannonly
namelist /exact/ test1,vl,vr,nv,gamma_law_eos
open(unit=11,file='modules/extriemanntest/problem.data',status='old',action='read')
read(unit=11,nml=exact)
close(11)
call initialize_global_parameters(gamma_law_eos,gamma_gas,'gamma_gas=')
mu=1.3d0
lscale=1.0d0
#if     puphase==1
    !generate the p-u phase
    open(unit=22,file='modules/extriemanntest/out/puphase.dat',status='replace',action='write')
    tl=test1(3)
    tr=test1(6)
    pl=test1(1)/(mu*amu)*kb*tl
    pr=test1(4)/(mu*amu)*kb*tr
    wl(1:2)=test1(1:2)
    wl(1)=wl(1)*(lscale**3)
    wl(2)=wl(2)/lscale
    wl(3)=pl*lscale
    wl(4)=gamma_gas
    wr(1:2)=test1(4:5)
    wr(1)=wr(1)*(lscale**3)
    wr(2)=wr(2)/lscale
    wr(3)=pr*lscale
    wr(4)=gamma_gas
    n_sample=1001
    dlogp=log10(wl(5)/8.525E-04)/(n_sample-1)
    print *,'puphase: ','pmax=',max(wl(5),wr(5)),'pmin=',min(wl(5),wr(5)),dlogp
    allocate(ustarl(n_sample),ustarr(n_sample),pstar_guess(n_sample))
    do i=1,n_sample
        pstar_guess(i)=8.525E-04*10**((i-1)*dlogp)
        ustarl(i)=exactrm_ustarl(pstar_guess(i),wl)
        ustarr(i)=exactrm_ustarr(pstar_guess(i),wr)
        write(22,fmt='(5ES15.4E2)') pstar_guess(i),ustarl(i),ustarr(i)
    end do
    deallocate(ustarl,ustarr,pstar_guess)
    close(22)
#else
    open(unit=23,file='modules/extriemanntest/out/extrm1.dat',status='replace',action='write')
    print *,gamma_gas
    !solve exact constant gamma gas riemann problem
    dv=(vr-vl)/(nv-1)
    boundl=0
    boundr=nv-1
    tolerance=1E-12
    do i=1,1
        tl=test1(3)
        tr=test1(6)
        pl=test1(1)/(mu*amu)*kb*tl
        pr=test1(4)/(mu*amu)*kb*tr
        do j=boundl,boundr
            t=1.0
            x=(vl+dv*j)*t
            wl(1:2)=test1(1:2)
            wl(1)=wl(1)*(lscale**3)
            wl(2)=wl(2)/lscale
            wl(3)=pl*lscale
            wl(4)=gamma_gas
            wr(1:2)=test1(4:5)
            wr(1)=wr(1)*(lscale**3)
            wr(2)=wr(2)/lscale
            wr(3)=pr*lscale
            wr(4)=gamma_gas
            call extrmsample(wl,wr,x,w,tolerance,t)
            temperature=w(3)/kb/(w(1)/(mu*amu))
            write(unit=23,fmt='(6ES18.6E2)') x/1e5,w(1),w(2)/1e5,w(3),temperature
            write(*,fmt='(5ES15.3E2)') x/1e5,w(1),w(2)/1e5,w(3),temperature
        end do
        close(23)
    end do
#endif
end program extrmtest
