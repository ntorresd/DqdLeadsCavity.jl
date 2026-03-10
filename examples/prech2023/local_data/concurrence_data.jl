N_range = 1000

tc_range = range(1e-5, 50. + 1e-5, N_range)
eV_range = range(0., 100., N_range)


begin
    local dqdObj = deepcopy(dqd_leads)
    C_matrix = Matrix{Float64}(undef, N_range, N_range)
    for (i, tc) in enumerate(tc_range)
        dqdObj.dqd.tc = tc
        print("tc = ", tc, "\n", N_range - i, " steps left\n")
        for (j, eV) in enumerate(eV_range)
            dqdObj.leads.Δμ = eV
            local ρss = steadystate(
                build_H_dqd_LR(dqdObj.dqd),
                build_L_ops_dqd_loc(dqdObj)
            )
            C_matrix[i, j] = get_concurrence_LR_num(dqdObj, ρss)
        end
    end
end