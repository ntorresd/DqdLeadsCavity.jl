using QuantumToolbox

N_points = 10000;
Ndot_range = logrange(1e-5, 1e3, N_points);

begin
    local dqdObj = deepcopy(dqd_leads_cavity);
    local ϵ_avg = dqdObj.dqd_leads.dqd.ϵ_avg;
    local μL, μR = get_chemical_potentials(dqdObj.dqd_leads);
    # analytical heat currents without the cavity
    J_ana_L_nd = [get_heat_current_gl(dqdObj.dqd_leads; left = true)]
    J_ana_R_nd = [get_heat_current_gl(dqdObj.dqd_leads; left = false)]
    # analytical heat current in the large-drive limit
    # TODO: find the right expression for this
    J_ana_L_ld = [get_heat_current_ld(dqdObj)]
    # numerical heat currents with the cavity
    JL_num_drive = [];
    JR_num_drive = [];
    for Ndot in Ndot_range
        dqdObj.cavity.Ndot = Ndot;
        local ρss = steadystate(
            build_H_z2022(dqdObj),
            build_L_ops_z2022(dqdObj)
        );
        # Operators
        local dL, dR = build_dqd_fermi_ops_LR(dqdObj);
        local N_op = dL' * dL + dR' * dR;
        local L_ops = build_L_ops_dqd_gl(dqdObj);
        local H_td = build_H_dqd_LR(dqdObj);
        local DLρ = D_sop(L_ops[1]) * ρss + D_sop(L_ops[2]) * ρss +
                    D_sop(L_ops[3]) * ρss + D_sop(L_ops[4]) * ρss;
        local DRρ = D_sop(L_ops[5]) * ρss + D_sop(L_ops[6]) * ρss +
                    D_sop(L_ops[7]) * ρss + D_sop(L_ops[8]) * ρss;

        push!(JL_num_drive, real(expect(H_td - μL * N_op, DLρ)));
        push!(JR_num_drive, real(expect(H_td - μR * N_op, DRρ)));
    end
end