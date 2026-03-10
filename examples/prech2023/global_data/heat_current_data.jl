using QuantumToolbox

N_points = 10000
tc_range = logrange(0.01, 1e3, N_points);

begin
	local dqdObj = deepcopy(dqd_leads);

    JL_ana_gl = [];
    J_ana_R_gl = [];
    JL_num_gl = [];
    J_num_R_gl = [];
    for tc in tc_range
        dqdObj.dqd.tc = tc;
        local L_ops = build_L_ops_dqd_gl(dqdObj);
        # steady state
		ρss = steadystate(
			build_H_dqd_ge(dqdObj.dqd),
			L_ops
		);
        # analytical heat currents
        push!(JL_ana_gl, get_heat_current_gl(dqdObj; left = true));
        push!(J_ana_R_gl, get_heat_current_gl(dqdObj; left = false));
        # numerical heat currents
        local H_td = build_H_dqd_ge(dqdObj.dqd);
        push!(
            JL_num_gl,
            get_heat_current_num(dqdObj, ρss, L_ops[1:4], H_td; left = true)
        );
        push!(
            J_num_R_gl,
            get_heat_current_num(dqdObj, ρss, L_ops[5:8], H_td; left = false)
        );
    end
end
