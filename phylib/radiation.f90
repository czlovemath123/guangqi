module radiation
#include <petsc/finclude/petscksp.h>
use phylib
use mathlib
use datastructure
use eos
use petsc_fld
!use petsc_fld_mg
use radiation_common_functions
implicit none

contains

subroutine radiation_transfer()
    real(8) :: t,dt,ddt,dt_remain
    integer :: i,ijk(3),ierr
    t=time_sys%t
    dt=time_sys%dt_hydro
    if (iradiation==1) then
    else if (iradiation==2) then
        !2d static short characteristics
    else if (iradiation==3) then
    else if (iradiation==4) then
        !petsc flux-limited diffusion
        if (refine_type=='static'.or.refine_type=='none') then
            call petsc_fld_solve(axb_operator,ierr)
        else if (refine_type=='adaptive'.or.refine_type=='mixed') then
            call petsc_fld_initialize(axb_operator,ierr)
            call petsc_fld_solve(axb_operator,ierr)
            call petsc_fld_recycle(axb_operator,ierr)
        end if
    else if (iradiation==5) then
    else
    end if
end subroutine radiation_transfer

subroutine initialize_radiation_environment()
    integer :: ierr
    if (iradiation==1) then
    else if (iradiation==2) then
    else if (iradiation==3) then
    else if (iradiation==4) then
        if (refine_type=='static'.or.refine_type=='none') then
            call generate_rosseland_planck_kap_table()
            call petsc_fld_initialize(axb_operator,ierr)
        else if (refine_type=='adaptive'.or.refine_type=='mixed') then
            call generate_rosseland_planck_kap_table()
        end if
    else if (iradiation==5) then
    end if
end subroutine initialize_radiation_environment

subroutine finalize_radiation()
    integer :: ierr
    if (iradiation==1) then
    else if (iradiation==2) then
    else if (iradiation==3) then
    else if (iradiation==4) then
        call petsc_fld_finalize(axb_operator,ierr)
    else if (iradiation==5) then
    end if
end subroutine finalize_radiation

end module radiation
