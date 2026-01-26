begin
    dqdObj = deepcopy(dqd_leads)
    results = map(tc_range) do tc
        dqdObj.dqd.tc = tc
        get_particle_current_neqgf(dqdObj)   # (J_L, J_R)
    end

    I_L_neqgf = first.(results)
    I_R_neqgf = last.(results)
end