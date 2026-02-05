using QuantumToolbox

N_points = 10000
tc_range = logrange(0.01, 1e3, N_points);

begin
    local dqdObj = deepcopy(dqd_leads)

    I_L_ana_gl = []
    I_R_ana_gl = []
    I_L_num_gl = []
    I_R_num_gl = []
    for tc in tc_range
        dqdObj.dqd.tc = tc
        # operators
        local L_ops = build_L_ops_dqd_gl(dqdObj)
        # steady state
		local ρss = steadystate(
			build_H_dqd_ge(dqdObj.dqd),
			L_ops
		)
        # analytical particle current
        push!(I_L_ana_gl, get_particle_current_gl(dqdObj; left = true))
        push!(I_R_ana_gl, get_particle_current_gl(dqdObj; left = false))
        # numerical particle currents
        push!(
            I_L_num_gl, 
            get_particle_current_num(dqdObj, ρss, L_ops[1:4])
        )
        push!(
            I_R_num_gl,
            get_particle_current_num(dqdObj, ρss, L_ops[5:8])
        )
    end
end
