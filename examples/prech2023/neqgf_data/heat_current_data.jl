N_points = 10000
tc_range = logrange(0.01, 1e3, N_points);

begin
    local tc_max = maximum(tc_range)
    local int_lims = (-10. * tc_max, 10. * tc_max);
    dqdObj = deepcopy(dqd_leads);

    JL_neqgf = [];
    J_R_neqgf = [];
    for tc in tc_range
        dqdObj.dqd.tc = tc
        push!(
            JL_neqgf,
            get_current_heat_neqgf(dqdObj; int_lims = int_lims, left = true)
        )
        push!(
            J_R_neqgf,
            get_current_heat_neqgf(dqdObj; int_lims = int_lims, left = false)
        )
    end
end