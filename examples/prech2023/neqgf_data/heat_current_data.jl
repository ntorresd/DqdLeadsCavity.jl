begin
    dqdObj = deepcopy(dqd_leads)
    results = map(tc_range) do tc
        dqdObj.dqd.tc = tc
        get_current_heat_neqgf(dqdObj)   # (J_L, J_R)
    end

    J_L_neqgf = first.(results)
    J_R_neqgf = last.(results)
end