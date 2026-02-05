using QuantumToolbox

N_points = 10000;
Ndot_range = logrange(1e-3, 1e3, N_points);

begin
    local dqdObj = deepcopy(dqd_leads_cavity);
    local dqdObj_ld = deepcopy(dqd_leads_cavity);
    local dqdObj_ld.cavity.Ndot = 5000;

    # Analytical large-drive particle current
    I_L_ld_ana = [get_particle_current_ld(dqdObj_ld.dqd_leads)];
    # Numerical large-drive particle current
    I_L_ld_num = [
        get_particle_current_num(
            dqdObj_ld,
            steadystate(
                build_H_z2022(dqdObj_ld),
                build_L_ops_z2022(dqdObj_ld)
            ),
            build_L_ops_dqd_gl(dqdObj_ld)[1:4]
        )
    ];
    # Numerical particle current
    I_L_num = [];
    for Ndot in Ndot_range
        dqdObj.cavity.Ndot = Ndot
        # steady state
        local ρss = steadystate(
            build_H_z2022(dqdObj),
            build_L_ops_z2022(dqdObj)
        );
        # numerical particle currents
        push!(
            I_L_num,
            get_particle_current_num(
                dqdObj,
                ρss,
                build_L_ops_dqd_gl(dqdObj)[1:4]
            )
        );
    end
end