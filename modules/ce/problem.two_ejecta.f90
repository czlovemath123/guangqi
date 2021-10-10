&rhd_quantities
    rho0=1d-15
    v0=0d0,0,0
    temp0=1000d0
    m_star=3.3e33
    t_outflow1=1d0       !outflow time span
    temp_outflow1=2.5d4    !outflow temperature
    mdot1=1d0           !solar mass per year
    v_esc_coefficient1=2.5d0
    t_outflow2=1d0       !outflow time span
    temp_outflow2=2.5d4    !outflow temperature
    mdot2=1d0           !solar mass per year
    v_esc_coefficient2=1.1d0
    dt_record=1d3
    petsc_rtol=1d-10
    record_length=13
    energy_conserve_formalism=.true.
    petsc_eos_rtol=2d-2
    resolve_rad_dt=.false.
    erad_slope_refine_thresh=0.02
    rtau=1d-2
/

