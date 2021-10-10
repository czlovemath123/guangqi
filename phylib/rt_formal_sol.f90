module rt_formal_sol
use datastructure
use mathlib
use phylib
implicit none

contains

function I_formal_sol(I0,source,chi,ds)
    !2nd order solution of radiation transfer equation
    !using Bezier polynomial interpolation
    real(8) :: I0,source(3),chi(3),ds(2),I_formal_sol,dtau,chi_control,phi(3)
    real(8) :: dtau_seg(2),source_control,smin,smax,chimin,chimax
    dtau_seg(1)=ds(1)/2d0*(chi(1)+chi(2))
    dtau_seg(2)=ds(2)/2d0*(chi(2)+chi(3))
    chi_control=f_control(chi,ds)
    chimin=min(chi(1),chi(2))
    chimax=max(chi(1),chi(2))
    if (chi(2)>chi(1).and.chi(2)>chi(3)) then
        !extremum at s0
        chi_control=chi(2)
    else if (chi_control<chimin.or.chi_control>chimax) then
        !overshoot
        chi_control=chi(1)
    end if
    dtau=ds(1)/3d0*(chi(1)+chi(2)+chi_control)
    source_control=f_control(source,dtau_seg)
    smin=min(source(1),source(2))
    smax=max(source(1),source(2))
    if (source(2)>source(1).and.source(2)>source(3)) then
        !extremum at s0
        phi(1)=extreme_coef1(dtau_seg(1))
        phi(2)=extreme_coef2(dtau_seg(1))
        phi(3)=0d0
    else if (source_control<smin.or.source_control>smax) then
        !overshoot
        phi(1)=overshoot_coef1(dtau_seg(1))
        phi(2)=overshoot_coef2(dtau_seg(1))
        phi(3)=0d0
    else
        phi(1)=bezier_poly_coef1(dtau_seg(1),dtau_seg(2))
        phi(2)=bezier_poly_coef2(dtau_seg(1),dtau_seg(2))
        phi(3)=bezier_poly_coef3(dtau_seg(1),dtau_seg(2))
    end if
    I_formal_sol=I0*exp(-dtau)+dot_product(phi,source)
end function I_formal_sol

function f_control(v,delta)
    !control point
    real(8) :: v(3),delta(2),f_control
    f_control=v(2)-delta(1)/2d0*(delta(2)/(delta(1)+delta(2))*(v(2)-v(1))/delta(1)  &
        +delta(1)/(delta(1)+delta(2))*(v(3)-v(2))/delta(2))
end function f_control

function bezier_aux_func1(t)
    real(8) :: t,bezier_aux_func1
    bezier_aux_func1=1d0-exp(-t)
end function bezier_aux_func1

function bezier_aux_func2(t)
    real(8) :: t,bezier_aux_func2
    bezier_aux_func2=t-bezier_aux_func1(t)
end function bezier_aux_func2

function bezier_aux_func3(t)
    real(8) :: t,bezier_aux_func3
    bezier_aux_func3=t**2-2d0*bezier_aux_func2(t)
end function bezier_aux_func3

function bezier_poly_coef1(dtau1,dtau2)
    !dtau1=delta tau upstream, dtau2=delta tau downstream
    real(8) :: dtau1,dtau2,bezier_poly_coef1
    bezier_poly_coef1=bezier_aux_func1(dtau1)+(bezier_aux_func3(dtau1)-(dtau2+2d0*dtau1))   &
        *bezier_aux_func2(dtau1)/dtau1/(dtau1+dtau2)
end function bezier_poly_coef1

function bezier_poly_coef2(dtau1,dtau2)
    real(8) :: dtau1,dtau2,bezier_poly_coef2
    bezier_poly_coef2=((dtau1+dtau2)*bezier_aux_func2(dtau1)-bezier_aux_func3(dtau1))/(dtau1*dtau2)
end function bezier_poly_coef2

function bezier_poly_coef3(dtau1,dtau2)
    real(8) :: dtau1,dtau2,bezier_poly_coef3
    bezier_poly_coef3=(bezier_aux_func3(dtau1)-dtau1*bezier_aux_func2(dtau1))/dtau2/(dtau1+dtau2)
end function bezier_poly_coef3

function extreme_coef1(dtau1)
    real(8) :: extreme_coef1,dtau1
    extreme_coef1=bezier_aux_func1(dtau1)+(bezier_aux_func3(dtau1)  &
        -2d0*dtau1*bezier_aux_func2(dtau1))/dtau1**2
end function extreme_coef1

function extreme_coef2(dtau1)
    real(8) :: extreme_coef2,dtau1
    extreme_coef2=(2d0*dtau1*bezier_aux_func2(dtau1)-bezier_aux_func3(dtau1))/dtau1**2
end function extreme_coef2

function overshoot_coef1(dtau1)
    real(8) :: overshoot_coef1,dtau1
    overshoot_coef1=bezier_aux_func1(dtau1)-bezier_aux_func3(dtau1)/dtau1**2
end function overshoot_coef1

function overshoot_coef2(dtau1)
    real(8) :: overshoot_coef2,dtau1
    overshoot_coef2=bezier_aux_func3(dtau1)/dtau1**2
end function overshoot_coef2

end module rt_formal_sol
