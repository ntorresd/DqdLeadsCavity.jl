export build_L_ops_local_LR, build_L_ops_semilocal_LR
export build_L_ops_local_ge

"""
Build Lindblad operators corresponding to the local LME evaluating
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
    cL, cR = build_dqd_ladder_ops_LR(dqd_leads.dqd)

    # Lindblad jump operators (lead-dot tunneling)
    L_ops = [
        sqrt(dqd_leads.ΓL * fL) * cL',
        sqrt(dqd_leads.ΓL * (1 - fL)) * cL,
        sqrt(dqd_leads.ΓR * fR) * cR',
        sqrt(dqd_leads.ΓR * (1 - fR)) * cR
    ]
    return L_ops
end

"""
Build Lindblad operators corresponding to the semi-local LME evaluating
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
    cL, cR = build_dqd_ladder_ops_LR(dqd_leads.dqd)

    # Lindblad jump operators (lead-dot tunneling)
    L_ops = [
        sqrt(ΓL * fL) * (1. - cR' * cR) * cL',
        sqrt(ΓL * (1. - fL)) * (1. - cR' * cR) * cL,
        sqrt(ΓL * fL_U) * cR' * cR * cL',
        sqrt(ΓL * (1. - fL_U)) * cR' * cR * cL,
        sqrt(ΓR * fR) * (1. - cL' * cL) * cR',
        sqrt(ΓR * (1. - fR)) * (1. - cL' * cL) * cR,
        sqrt(ΓR * fR_U) * cL' * cL * cR',
        sqrt(ΓR * (1. - fR_U)) * cL' * cL * cR,
    ]
    return L_ops
end

"""
Build Lindblad operators corresponding to the local LME evaluating
the eigen-energies of the DQD (g-e)
"""
function build_L_ops_local_ge(dqd_leads::DqdLeads)
    # Parameters
    Ω = get_Ω(dqd_leads.dqd); θ = get_θ(dqd_leads.dqd)
    ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
    Δμ = dqd_leads.leads.Δμ

    # Fermi functions
    fLg = fermi(ϵg, Δμ/2., TL)
    fRg = fermi(ϵg, - Δμ/2., TR)

    fLe = fermi(ϵe, Δμ/2., TL)
    fRe = fermi(ϵe, - Δμ/2., TR)

    # Jump rates by process
    Γg0L = ΓL * fLg * cos(θ/2)^2
    Γg0R = ΓR * fRg * sin(θ/2)^2
    Γe0L = ΓL * fLe * sin(θ/2)^2
    Γe0R = ΓR * fRe * cos(θ/2)^2
    Γ0gL = ΓL * (1. - fLg) * cos(θ/2)^2
    Γ0gR = ΓR * (1. - fRg) * sin(θ/2)^2
    Γ0eL = ΓL * (1. - fLe) * sin(θ/2)^2
    Γ0eR = ΓR * (1. - fRe) * cos(θ/2)^2
    
    # Creation and annihilation operators 
    cg, ce = build_dqd_ladder_ops_ge(dqd_leads.dqd)

    # Lindblad operators
    L_ops = [
        sqrt(Γg0L) * cg', sqrt(Γ0gL) * cg, sqrt(Γe0L) * ce',sqrt(Γ0eL) * ce,
        sqrt(Γg0R) * cg', sqrt(Γ0gR) * cg, sqrt(Γe0R) * ce', sqrt(Γ0eR) * ce
    ]
    return L_ops
end