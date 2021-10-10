module source_control
use datastructure
use mathlib
use eos
use phylib
use problem
use geometricsource
use gravity
use cooling
implicit none

contains

subroutine source_apply(blk)
    !currently, only cooling, geometry and gravity are calculated by source_apply
    type(blockdef), pointer :: blk
    character(len=128) :: alert
    integer :: ijk(3),i,j,k
    real(8) :: t,dt,dedt,de,t_source,t_source_final,lum_atm,rthetaphi(3),ml_rate,rho
    real(8), dimension(5) :: u_temp,source,dsource,source_out,dsource_stiff,dsourcedt
    !need to calculate additional quantities when specific physics is use, e.g. divv and opacity
    dt=time_sys%dt_hydro
    t=time_sys%t
    if (nd==1) then
    else if (nd==2) then
        if (igravity/=0) then
            call gravity_source(blk)
            call convert_u_to_source_block(blk)
            blk%source=blk%source+blk%s_grav*dt
            call convert_source_to_u_block(blk)
            call convert_u_to_w_block(blk)
        end if
    end if
end subroutine source_apply

end module source_control
