using QuantumToolbox

begin
	dqdObj = deepcopy(dqd_leads)
	ρss_list = []
	for tc in tc_range
		dqdObj.dqd.tc = tc
		ρss = steadystate(
			build_H_dqd_ge(dqdObj.dqd),
			build_L_ops_dqd_gl(dqdObj)
		)
		push!(ρss_list, ρss)
	end
end

begin
    J_ana_L_gl = []
    J_ana_R_gl = []

    J_num_L_gl = []
    J_num_R_gl = []

    local i = 1
    for tc in tc_range
        dqdObj.dqd.tc = tc
        local μL, μR = get_chemical_potentials(dqdObj.leads)
        local H_dqd = build_H_dqd_LR(dqdObj.dqd)
        local dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
        local N_op = dL' * dL + dR' * dR
        local L_ops = build_L_ops_dqd_gl(dqdObj)
        # analytical heat currents
        local J_L, J_R = get_heat_current_gl(dqdObj)
        push!(J_ana_L_gl, J_L)
        push!(J_ana_R_gl, J_R)
        # numerical heat currents
        local DLρ = D_sop(L_ops[1], ρss_list[i]) + D_sop(L_ops[2], ρss_list[i]) +
              D_sop(L_ops[3], ρss_list[i]) + D_sop(L_ops[4], ρss_list[i])
        local DRρ = D_sop(L_ops[5], ρss_list[i]) + D_sop(L_ops[6], ρss_list[i]) +
              D_sop(L_ops[7], ρss_list[i]) + D_sop(L_ops[8], ρss_list[i])
        push!(J_num_L_gl, real(expect(H_dqd - μL * N_op, DLρ)))
        push!(J_num_R_gl, real(expect(H_dqd - μR * N_op, DRρ)))
        i += 1
    end
end
