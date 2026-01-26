using QuantumToolbox

begin
	dqdObj = deepcopy(dqd_leads)
	ρss_list = []
	for tc in tc_range
		dqdObj.dqd.tc = tc
		ρss = steadystate(
			build_H_dqd_LR(dqdObj.dqd),
			build_L_ops_dqd_loc(dqdObj)
		)
		push!(ρss_list, ρss)
	end
end

begin
    J_ana_L_loc = []
    J_ana_R_loc = []

    J_num_L_loc = []
    J_num_R_loc = []

    local i = 1
    for tc in tc_range
        dqdObj.dqd.tc = tc
        local ϵ_avg = dqdObj.dqd.ϵ_avg
        local μL, μR = get_chemical_potentials(dqdObj.leads)
        local dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
        local N_op = dL' * dL + dR' * dR
        local L_ops = build_L_ops_dqd_loc(dqdObj)
        # analytical heat currents
        local J_L, J_R = get_heat_current_loc(dqdObj)
        push!(J_ana_L_loc, J_L)
        push!(J_ana_R_loc, J_R)
        # numerical heat currents
        local H_td = ϵ_avg * (dL' * dL + dR' * dR)
        local DLρ = D_sop(L_ops[1], ρss_list[i]) + D_sop(L_ops[2], ρss_list[i])
        local DRρ = D_sop(L_ops[3], ρss_list[i]) + D_sop(L_ops[4], ρss_list[i])
        push!(J_num_L_loc, real(expect(H_td - μL * N_op, DLρ)))
        push!(J_num_R_loc, real(expect(H_td - μR * N_op, DRρ)))
        i += 1
    end
end
