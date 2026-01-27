begin
    local tc_max = maximum(tc_range)
    dqdObj = deepcopy(dqd_leads)
    results = map(tc_range) do tc
        dqdObj.dqd.tc = tc
        get_current_heat_neqgf(dqdObj; int_lims = (-10. * tc_max, 10. * tc_max))   # (J_L, J_R)
    end

    J_L_neqgf = first.(results)
    J_R_neqgf = last.(results)
end