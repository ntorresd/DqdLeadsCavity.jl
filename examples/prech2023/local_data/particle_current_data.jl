using QuantumToolbox

N_points = 10000
tc_range = logrange(0.01, 1e3, N_points);

begin
	local dqdObj = deepcopy(dqd_leads)

    I_L_ana_loc = []
    I_R_ana_loc = []
    I_L_num_loc = []
    I_R_num_loc = []
    for tc in tc_range
        dqdObj.dqd.tc = tc
        # operators
        local L_ops = build_L_ops_dqd_loc(dqdObj)
        # steady state
		local ρss = steadystate(
			build_H_dqd_LR(dqdObj.dqd),
			L_ops
		)
        # analytical particle current
        push!(I_L_ana_loc, get_particle_current_loc(dqdObj; left = true))
        push!(I_R_ana_loc, get_particle_current_loc(dqdObj; left = false))
        # numerical particle currents
        push!(
            I_L_num_loc,
            get_particle_current_num(dqdObj, ρss, L_ops[1:2])
        )
        push!(
            I_R_num_loc,
            get_particle_current_num(dqdObj, ρss, L_ops[3:4])
        )
    end
end
