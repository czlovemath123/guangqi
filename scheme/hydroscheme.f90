module hydroscheme
#if     ischeme==0
    !has no hydrodynamics, fluid does not evolve
    use godunov
#elif   ischeme==1
    use godunov
#elif   ischeme==2
    use muscl
#endif

implicit none
contains

subroutine hydro_step()
    call communicate_hydro()
    call applyboundconds()
#if     ischeme==0
    time_sys%dt_hydro=time_sys%dt_hydro*1.01
    if (time_sys%dt_hydro+time_sys%t>time_sys%t_final) then
        time_sys%dt_hydro=time_sys%t_final-time_sys%t
    end if
#elif   ischeme==1
    call godunov_hydro_split()
#elif   ischeme==2
    call vanleer_hydro_unsplit()
#endif
end subroutine hydro_step

end module hydroscheme
