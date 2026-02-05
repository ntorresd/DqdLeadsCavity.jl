begin
    local tc_max = maximum(tc_range)
    local int_lims = (-10. * tc_max, 10. * tc_max)
    dqdObj = deepcopy(dqd_leads)

    I_L_neqgf = []
    I_R_neqgf = []
    for tc in tc_range
        dqdObj.dqd.tc = tc
        push!(
            I_L_neqgf,
            get_particle_current_neqgf(dqdObj; left = true, int_lims = int_lims)
        )
        push!(
            I_R_neqgf,
            get_particle_current_neqgf(dqdObj; left = false, int_lims = int_lims)
        )
    end
end