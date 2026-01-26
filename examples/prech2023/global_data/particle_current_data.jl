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
    I_L_ana_gl = []
    I_R_ana_gl = []

    I_L_num_gl = []
    I_R_num_gl = []

    local i = 1
    for tc in tc_range
        dqdObj.dqd.tc = tc
        local μL, μR = get_chemical_potentials(dqdObj.leads)
        local H_dqd = build_H_dqd_LR(dqdObj.dqd)
        local dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
        local N_op = dL' * dL + dR' * dR
        local L_ops = build_L_ops_dqd_gl(dqdObj)
        # analytical particle current
        local I_L, I_R = get_particle_current_gl(dqdObj)
        push!(I_L_ana_gl, I_L)
        push!(I_R_ana_gl, I_R)
        # numerical heat currents
        local DLρ = D_sop(L_ops[1], ρss_list[i]) + D_sop(L_ops[2], ρss_list[i]) +
              D_sop(L_ops[3], ρss_list[i]) + D_sop(L_ops[4], ρss_list[i])
        local DRρ = D_sop(L_ops[5], ρss_list[i]) + D_sop(L_ops[6], ρss_list[i]) +
              D_sop(L_ops[7], ρss_list[i]) + D_sop(L_ops[8], ρss_list[i])
        push!(I_L_num_gl, real(expect(N_op, DLρ)))
        push!(I_R_num_gl, real(expect(N_op, DRρ)))
        i += 1
    end
end
