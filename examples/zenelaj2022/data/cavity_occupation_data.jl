using QuantumToolbox

N_points = 10000;
Ndot_range = range(0.0, 25.0, N_points);

begin
    local dqdObj = deepcopy(dqd_leads_cavity);

    n_cav_ss_num = [];
    n_cav_ss_ana = [];
    for Ndot in Ndot_range
        dqdObj.cavity.Ndot = Ndot;
        local ρss = steadystate(
            build_H_z2022(dqdObj),
            build_L_ops_z2022(dqdObj)
        );
        local a = build_cav_a_op(dqdObj; disp = true);

        push!(n_cav_ss_num, real(expect(a' * a, ρss)));
        push!(n_cav_ss_ana, get_n_cav_ss(dqdObj.cavity));
    end
end