export build_H_dqd_LR, build_H_dqd_ge

@doc raw"""
Build DQD Hamiltonian in the L-R basis
"""
function build_H_dqd_LR(dqd::Dqd)
    # Parameters
    ϵL, ϵR = get_onsite_energies(dqd)
    # Operators
    cL, cR = build_dqd_ladder_ops_LR(dqd)
    nL, nR, nD = build_dqd_number_ops_LR(dqd)


    # Hamiltonian: on-site + tunneling
    H = ϵL * cL' * cL + ϵR * cR' * cR + dqd.tc * (cL' * cR + cR' * cL)

    # Add non-linear term for finite Coulomb interaction
    H = dqd.blockade ? H : H + dqd.U * cL' * cL * cR' * cR

    return H
end


@doc raw"""
Build σz in the g-e parametrization
"""
function build_σz_ge(dqd::Dqd)
    ket_0, ket_g, ket_e = build_dqd_basis_ge(dqd)
    σz = ket_e * ket_e' - ket_g * ket_g'
    return σz
end

@doc raw"""
Build DQD Hamiltonian in the g-e parametrization
"""
function build_H_dqd_ge(dqd::Dqd)
    cg, ce = build_dqd_ladder_ops_ge(dqd)
    ϵg, ϵe = get_eigen_energies(dqd)
    H = (ϵg * cg' * cg) + (ϵe * ce' * ce)

    # Add non-linear term for finite Coulomb interaction
    H = dqd.blockade ? H : H + (dqd.U * cg' * cg * ce' * ce)
    return H
end
function build_H_dqd_ge(dqd::Dqd, cavity::Cavity)
    dqd.blockade || throw(
        ArgumentError(
            "DQD g-e basis representation is not supported for dqd.blockade=$(dqd.blockade)"
        ),
    )
    # Parameters
    dim = get_dim(dqd)
    Δd = get_Δd(dqd, cavity)

    # Operators
    ## Identity with zero in the first entry (|0><0|)
    id_ge = id_no_vacuum(dim)
    σz_ge = build_σz_ge(dqd)

    H_dqd = Δd * σz_ge / 2. + dqd.ϵ_avg * id_ge
    return H_dqd
end