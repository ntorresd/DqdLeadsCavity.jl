using QuantumToolbox

N_points = 10000
tc_range = logrange(0.01, 1e3, N_points);

begin
	dqdObj = deepcopy(dqd_leads)

    J_ana_L_loc = []
    J_ana_R_loc = []
    J_num_L_loc = []
    J_num_R_loc = []
    for tc in tc_range
        # parameters
        dqdObj.dqd.tc = tc
        ϵ_avg = dqdObj.dqd.ϵ_avg
		# operators
        local L_ops = build_L_ops_dqd_loc(dqdObj);
        local H_dqd = build_H_dqd_LR(dqdObj.dqd);
        local ρss = steadystate(H_dqd, L_ops);
        # bookeeping Hamiltonian
        local dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
        local H_td = ϵ_avg * (dL' * dL + dR' * dR)
        # analytical heat currents
        push!(J_ana_L_loc, get_heat_current_loc(dqdObj; left = true))
        push!(J_ana_R_loc, get_heat_current_loc(dqdObj; left = false))
        # numerical heat currents
        push!(
            J_num_L_loc,
            get_heat_current_num(dqdObj, ρss, L_ops[1:2], H_td; left = true)
        )
        push!(
            J_num_R_loc,
            get_heat_current_num(dqdObj, ρss, L_ops[3:4], H_td; left = false)
        )
    end
end
