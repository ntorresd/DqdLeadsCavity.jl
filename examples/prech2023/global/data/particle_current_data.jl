using QuantumToolbox

begin
	tc_list = logrange(1e-3, 1e1, 1000);
	dqdObj = deepcopy(dqd_leads)
	ρss_list = []
	for tc in tc_list
		dqdObj.dqd.tc = tc
		ρss = steadystate(
			build_H_dqd_ge(dqdObj.dqd),
			build_L_ops_dqd_gl(dqdObj)
		)
		push!(ρss_list, ρss)
	end
end

begin
    I_ana = []
    I_num = []

    local i = 1
    for tc in tc_list
        dqdObj.dqd.tc = tc
        μL, μR = get_chemical_potentials(dqdObj.leads)
        H_dqd = build_H_dqd_LR(dqdObj.dqd)
        dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
        N_op = dL' * dL + dR' * dR
        L_ops = build_L_ops_dqd_gl(dqdObj)
        # analytical particle current
        push!(I_ana, get_particle_current_gl(dqdObj))
        # numerical heat currents
        DLρ = D_sop(L_ops[1], ρss_list[i]) + D_sop(L_ops[2], ρss_list[i]) +
              D_sop(L_ops[3], ρss_list[i]) + D_sop(L_ops[4], ρss_list[i])
        DRρ = D_sop(L_ops[5], ρss_list[i]) + D_sop(L_ops[6], ρss_list[i]) +
              D_sop(L_ops[7], ρss_list[i]) + D_sop(L_ops[8], ρss_list[i])
        push!(I_num, real(expect(N_op, DLρ)))
        i += 1
    end
end
