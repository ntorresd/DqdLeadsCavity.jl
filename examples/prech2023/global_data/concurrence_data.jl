N_range = 500

Δϵ_logrange = logrange(10, 1e4, N_range)
eV_logrange = logrange(1e1, 1e4, N_range)

begin
    local dqdObj = deepcopy(dqd_leads)
    C_list_Δϵ = Vector{Float64}(undef, length(Δϵ_logrange))
    for (i, Δϵ) in enumerate(Δϵ_logrange)
        dqdObj.dqd.Δϵ = Δϵ
        C_list_Δϵ[i] = get_concurrence_gl_ana(dqdObj)
    end
    
    # Concurrence as a function of eV for the detuning that gives
    # max concurrence for dqd_leads.leads.Δμ
    # (C_max, idx) = findmax(C_list_Δϵ)
    # dqdObj.dqd.Δϵ = Δϵ_logrange[idx]
    dqdObj = deepcopy(dqd_leads)
    C_list_eV = Vector{Float64}(undef, length(eV_logrange))
    for (j, eV) in enumerate(eV_logrange)
        dqdObj.leads.Δμ = eV
        C_list_eV[j] = get_concurrence_gl_ana(dqdObj)
    end
end

begin
    local dqdObj = deepcopy(dqd_leads)

    C_matrix_ana = Matrix{Float64}(undef, length(Δϵ_logrange), length(eV_logrange))
    C_matrix_num = Matrix{Float64}(undef, length(Δϵ_logrange), length(eV_logrange))
    pD_matrix_num = Matrix{Float64}(undef, length(Δϵ_logrange), length(eV_logrange))
    α_matrix_num = Matrix{Float64}(undef, length(Δϵ_logrange), length(eV_logrange))

    for (i, Δϵ) in enumerate(Δϵ_logrange)
        dqdObj.dqd.Δϵ = Δϵ
        print("Δϵ = ", Δϵ, "\n", N_range - i, " steps left\n")
        for (j, eV) in enumerate(eV_logrange)
            dqdObj.leads.Δμ = eV
            # analytical concurrence
            C_matrix_ana[i, j] = get_concurrence_gl_ana(dqdObj)

            # numerical concurrence
            local ρss =  steadystate(
                build_H_dqd_ge(dqdObj.dqd),
                build_L_ops_dqd_gl(dqdObj)
            )
            local ket_0, ket_L, ket_R, ket_D = build_dqd_basis_LR(dqdObj.dqd)
            C_matrix_num[i, j] = get_concurrence_LR_num(dqdObj, ρss)
            pD_matrix_num[i, j] = abs(ket_D' * ρss * ket_D)
            α_matrix_num[i, j] = get_coherence_LR_num(dqdObj, ρss)
        end
    end
end
