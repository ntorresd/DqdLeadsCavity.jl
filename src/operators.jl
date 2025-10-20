export build_H_dqd_LR, build_H_dqd_ge

"""
Build DQD Hamiltonian in the L-R basis
"""
function build_H_dqd_LR(dqd::Dqd)
    # Parameters
    ϵL, ϵR = get_onsite_energies(dqd)
    # Operators
    cL, cR = build_dqd_ladder_ops_LR(dqd)
    nL, nR, nD = build_dqd_number_ops_LR(dqd)


    # Hamiltonian: on-site + tunneling
    H = ϵL * nL + ϵR * nR + dqd.tc * (cL' * cR + cR' * cL)

    # Add non-linear term for finite Coulomb interaction
    H = dqd.blockade ? H : H + dqd.U * nL * nR

    return H
end

"""
Build DQD Hamiltonian in the g-e basis
"""
function build_H_dqd_ge(dqd::Dqd)
    dqd.blockade || throw(
        ArgumentError(
            "DQD g-e basis representation is not supported for dqd.blockade=$(dqd.blockade)"
        ),
    )
    # Parameters
    dim = get_dim(dqd)
    Ω = get_Ω(dqd)

    # Operators
    cg, ce = build_dqd_ladder_ops_ge(dqd)
    σz = ce' * ce - cg' * cg
    ## Identity with zero in the first entry (|0><0|)
    Id_ge = id_no_vacuum(dim)

    H_dqd = Ω * σz / 2. + dqd.ϵ_avg * Id_ge
    return H_dqd
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
    cg, ce = build_dqd_ladder_ops_ge(dqd)
    σz = ce' * ce - cg' * cg
    ## Identity with zero in the first entry (|0><0|)
    Id_ge = Id_ge = id_no_vacuum(dim)

    H_dqd = Δd * σz / 2. + dqd.ϵ_avg * Id_ge
    return H_dqd
end