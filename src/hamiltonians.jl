export build_H_dqd_LR, build_H_dqd_ge
export build_H_cav_drive_coh, build_H_JC

@doc raw"""
DQD Hamiltonian in the L-R basis
"""
function build_H_dqd_LR(dqd::Dqd)
    # Parameters
    ϵL, ϵR = get_onsite_energies(dqd)

    # Operators
    dL, dR = dqd.blockade ? build_dqd_vladder_ops_LR(dqd) : (fdestroy(2, 1), fdestroy(2, 2))

    # Hamiltonian: on-site + tunneling
    H = ϵL * dL' * dL + ϵR * dR' * dR + dqd.tc * (dL' * dR + dR' * dL)

    # Add non-linear term for finite Coulomb interaction
    H = dqd.blockade ? H : H + dqd.U * dL' * dL * dR' * dR

    return H
end


@doc raw"""
σz in the g-e parametrization
"""
function build_σz_ge(dqd::Dqd)
    ket_0, ket_g, ket_e = build_dqd_basis_ge(dqd)
    σz = ket_e * ket_e' - ket_g * ket_g'
    return σz
end

@doc raw"""
DQD Hamiltonian in the g-e parametrization
"""
function build_H_dqd_ge(dqd::Dqd)
    dg, de = dqd.blockade ? build_dqd_vladder_ops_ge(dqd) : build_dqd_fermi_ops_ge(dqd)
    ϵg, ϵe = get_eigen_energies(dqd)
    H = (ϵg * dg' * dg) + (ϵe * de' * de)

    # Add non-linear term for finite Coulomb interaction
    H = dqd.blockade ? H : H + (dqd.U * dg' * dg * de' * de)
    return H
end
function build_H_dqd_ge(dqd_leads_cavity::DqdLeadsCavityObj)
    dqd = dqd_leads_cavity.dqd_leads.dqd
    cavity = dqd_leads_cavity.cavity

    dqd.blockade || throw(
        ArgumentError(
            "DQD g-e basis representation is not supported for dqd.blockade=$(dqd.blockade)"
        ),
    )
    # Parameters
    Δd = get_Δd(dqd_leads_cavity)

    # Operators
    ## Identity with zero in the first entry (|0><0|)
    id_ge = tensor(build_id_dqd(dqd), build_id_cavity(cavity))
    σz_ge = build_dqd_σz_op(dqd_leads_cavity)

    H_dqd = Δd * σz_ge / 2. + dqd.ϵ_avg * id_ge
    return H_dqd
end

@doc raw"""
Hamiltonian for the coherent cavity drive
"""
function build_H_cav_drive_coh(cavity::Cavity; disp::Bool = true)
    # Parameters
    Δc = get_Δc(cavity)

    # Operators
    a = build_cav_a_op(cavity, disp = disp)

    # Hamiltonian
    H_drive = Δc * a' * a
    return H_drive
end
function build_H_cav_drive_coh(dqd_leads_cavity::DqdLeadsCavityObj; disp::Bool = true)
    # Parameters
    κin = dqd_leads_cavity.cavity.κin
    Ndot = dqd_leads_cavity.cavity.Ndot
    Δc = get_Δc(dqd_leads_cavity.cavity)

    # Operators
    a = build_cav_a_op(dqd_leads_cavity; disp = disp)

    # Hamiltonian
    H_drive = Δc * a' * a + sqrt(κin * Ndot) * (a' + a)
    return H_drive
end

@doc raw"""
Hamiltonian for the DQD-cavity interaction in the JC model
"""
function build_H_JC(dqd_leads_cavity::DqdLeadsCavityObj; disp::Bool = true)
    # Parameters
    gJC = dqd_leads_cavity.cavity.gJC

    # Operators
    a = build_cav_a_op(dqd_leads_cavity; disp = disp)
    σm, σp = build_dqd_ladder_ops_ge(dqd_leads_cavity)

    # Hamiltonian
    H_JC = gJC * (a' * σm + a * σp)
    return H_JC
end