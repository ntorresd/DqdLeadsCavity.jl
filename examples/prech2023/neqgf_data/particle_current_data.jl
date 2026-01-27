begin
    local tc_max = maximum(tc_range)
    dqdObj = deepcopy(dqd_leads)
    results = map(tc_range) do tc
        dqdObj.dqd.tc = tc
        get_particle_current_neqgf(dqdObj; int_lims = (-10. * tc_max, 10. * tc_max))   # (I_L, I_R)
    end

    I_L_neqgf = first.(results)
    I_R_neqgf = last.(results)
end