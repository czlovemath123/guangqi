module limiters
use mathlib

implicit none

contains

function slope_limiter(r)
    real(8) :: slope_limiter,omega,r,xi_l,xi_r,rnew
    if (r<=zero) then
        slope_limiter=zero
    else
        omega=0d0
        slope_limiter=min(2/(1+r),2*r/(1+r))
    end if
end function slope_limiter

function find_the_slope(diff_left,diff_right,diff_central)
    real(8) :: diff_left,diff_right,diff_central,find_the_slope,r
#if     isolver==1
    if (diff_central/=zero) then
        if (diff_left*diff_right<=0d0) then
            find_the_slope=0d0
        else
            find_the_slope=abs(diff_central)/diff_central*min(abs(diff_left),abs(diff_right))
            if (igeometry==1) then
                find_the_slope=find_the_slope*1d0
            else if (igeometry==2) then
                find_the_slope=find_the_slope*0.75d0
            end if
        end if
    else
        find_the_slope=0d0
    end if
#elif   isolver==2
    if (diff_central/=zero) then
        if (diff_left*diff_right<=0) then
            find_the_slope=0d0
        else
            find_the_slope=abs(diff_central)/diff_central*0.75d0*min(abs(diff_left),abs(diff_right))
        end if
    else
        find_the_slope=0d0
    end if
#endif
end function find_the_slope

function minmod(diff_left,diff_right)
    real(8) :: diff_left,diff_right,minmod
    if (diff_left*diff_right<=0) then
        minmod=0d0
    else
        minmod=abs(diff_left)/diff_left*min(abs(diff_left),abs(diff_right))
    end if
end function minmod

end module limiters
