export build_L_ops_local_LR, build_L_ops_semilocal_LR
export build_L_ops_thcg_int


@doc raw"""
Lindblad operators corresponding to the local LME evaluating
the fermi distributions at the corresponding on-site energies.
"""
function build_L_ops_local_LR(dqd_leads::DqdLeads)
    # Parameters
    ϵL, ϵR = get_onsite_energies(dqd_leads.dqd)
    μL, μR = get_chemical_potentials(dqd_leads.leads)

    # Fermi functions
    fL = fermi(ϵL, μL, dqd_leads.leads.TL)
    fR = fermi(ϵR, μR, dqd_leads.leads.TR)

    # Creation and annihilation operators 
    dL, dR = build_dqd_fermi_ops_LR(dqd_leads.dqd)

    # Lindblad jump operators (lead-dot tunneling)
    L_ops = [
        sqrt(dqd_leads.ΓL * fL) * dL',
        sqrt(dqd_leads.ΓL * (1 - fL)) * dL,
        sqrt(dqd_leads.ΓR * fR) * dR',
        sqrt(dqd_leads.ΓR * (1 - fR)) * dR
    ]
    return L_ops
end

@doc raw"""
Lindblad operators corresponding to the semi-local LME evaluating
the fermi distributions at the corresponding on-site energies.
"""
function build_L_ops_semilocal_LR(dqd_leads::DqdLeads)
    !dqd_leads.dqd.blockade || throw(
        ArgumentError(
            "Semi-local LME requires dqd.blockade=false and finite dqd.U"
        ),
    )
    # Parameters
    ϵL, ϵR = get_onsite_energies(dqd_leads.dqd)
    μL, μR = get_chemical_potentials(dqd_leads.leads)
    ΓL = dqd_leads.ΓL; ΓR = dqd_leads.ΓR
    U = dqd_leads.dqd.U

    # Fermi Functions
    fL = fermi(ϵL, μL, TL)
    fR = fermi(ϵR, μR, TR)
    fL_U = fermi(ϵL + U, μL, TL)
    fR_U = fermi(ϵR + U, μR, TR)
    
    # Creation and annihilation operators 
    dL, dR = build_dqd_fermi_ops_LR(dqd_leads.dqd)

    # Lindblad jump operators (lead-dot tunneling),
    L_ops = [
        sqrt(ΓL * fL) * (1. - dR' * dR) * dL',
        sqrt(ΓL * (1. - fL)) * (1. - dR' * dR) * dL,
        sqrt(ΓL * fL_U) * dR' * dR * dL',
        sqrt(ΓL * (1. - fL_U)) * dR' * dR * dL,
        sqrt(ΓR * fR) * (1. - dL' * dL) * dR',
        sqrt(ΓR * (1. - fR)) * (1. - dL' * dL) * dR,
        sqrt(ΓR * fR_U) * dL' * dL * dR',
        sqrt(ΓR * (1. - fR_U)) * dL' * dL * dR,
    ]
    return L_ops
end

@doc raw"""
Lindblad operators for the interacting DQD for the THC global case 
according to eq. (118) [Potts2021]. This only works for ϵL = ϵR.
"""
function build_L_ops_thcg_int(dqd_leads::DqdLeads)
    (; ϵL, ϵR, μL, μR, ΓL, ΓR, TL, TR) = context_LR(dqd_leads)
    ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
    U = dqd_leads.dqd.U

    fLg, fLe, fRg, fRe = get_fermi_ge(dqd_leads)
    fLgU, fLeU, fRgU, fReU = fermi(ϵg + U, μL, TL), fermi(ϵe + U, μL, TL), fermi(ϵg + U, μR, TR), fermi(ϵe + U, μR, TR)

    cg, ce = build_dqd_ladder_ops_ge(dqd_leads.dqd)
    L_ops = [
        # Left lead jump operators
        ## ground state
        sqrt(ΓL * fLg / 2.) * (1. - ce' * ce) * cg,
        sqrt(ΓL * (1. - fLg) / 2.) * (1. - ce' * ce) * cg',
        sqrt(ΓL * fLgU / 2.) * (ce' * ce * cg'),
        sqrt(ΓL * (1. - fLgU) / 2.) * (ce' * ce * cg),
        ## excited state
        sqrt(ΓL * fLe / 2.) * (1. - cg' * cg) * ce,
        sqrt(ΓL * (1. - fLe) / 2.) * (1. - cg' * cg) * ce',
        sqrt(ΓL * fLeU / 2.) * (cg' * cg * ce'),
        sqrt(ΓL * (1. - fLeU) / 2.) * (cg' * cg * ce),
        # Right lead jump operators
        ## ground state
        sqrt(ΓR * fRg / 2.) * (1. - ce' * ce) * cg,
        sqrt(ΓR * (1. - fRg) / 2.) * (1. - ce' * ce) * cg',
        sqrt(ΓR * fRgU / 2.) * (ce' * ce * cg'),
        sqrt(ΓR * (1. - fRgU) / 2.) * (ce' * ce * cg),
        ## excited state
        sqrt(ΓR * fRe / 2.) * (1. - cg' * cg) * ce,
        sqrt(ΓR * (1. - fRe) / 2.) * (1. - cg' * cg) * ce',
        sqrt(ΓR * fReU / 2.) * (cg' * cg * ce'),
        sqrt(ΓR * (1. - fReU) / 2.) * (cg' * cg * ce)
    ]

    return(L_ops)
end
