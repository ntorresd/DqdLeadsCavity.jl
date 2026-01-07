export current_heat_avg_thcg_ana
export build_L_ops_thcg

@doc raw"""
Calculates analytical steady-state heat current through the non-interacting DQD
for the thermodynamically consistent global LME according to eq. (B.8) [Potts et. al. 2021]

# Returns
- Steady-state average of the current operator

See also: [`build_dqd_basis`](@ref)
"""
function current_heat_avg_thcg_ana(dqd_leads::DqdLeads; left::Bool = true)
    # Parameters g-e
    (; ϵg, ϵe, ΓLg, ΓLe, ΓRg, ΓRe, fLg, fLe, fRg, fRe) = context_ge(dqd_leads)
    μL, μR = get_chemical_potentials(dqd_leads)

    Δfg = left ? fLg - fRg : fRg - fLg
    Δfe = left ? fLe - fRe : fRe - fLe
    μ = left ? μL : μR

    γg = ΓLg * ΓRg / (ΓLg + ΓRg)
    γe = ΓLe * ΓRe / (ΓLe + ΓRe)

    Jg = γg * (ϵg - μ) * Δfg
    Je = γe * (ϵe - μ) * Δfe

    Jα = Jg + Je
    return Jα
end

@doc raw"""
Build Lindblad operators corresponding to the local LME evaluating
the eigen-energies of the DQD (g-e)
"""
function build_L_ops_thcg(dqd_leads::DqdLeads)
    # Parameters g-e
    (; ϵg, ϵe, ΓLg, ΓLe, ΓRg, ΓRe, fLg, fLe, fRg, fRe) = context_ge(dqd_leads)
    μL, μR = get_chemical_potentials(dqd_leads)
    TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR
    # Fermi functions
    fLg = fermi(ϵg, μL, TL)
    fRg = fermi(ϵg, μR, TR)

    fLe = fermi(ϵe, μL, TL)
    fRe = fermi(ϵe, - μR, TR)
    
    # Creation and annihilation operators 
    cg, ce = build_dqd_ladder_ops_ge(dqd_leads.dqd)

    # Lindblad operators
    L_ops = [
        sqrt(ΓLg * fLg) * cg' + sqrt(ΓLg * (1. - fLg)) * cg,
        sqrt(ΓLe * fLe) * ce' + sqrt(ΓLe * (1. - fLe)) * ce,
        sqrt(ΓRg * fRg) * cg' + sqrt(ΓRg * (1. - fRg)) * cg,
        sqrt(ΓRe * fRe) * ce' + sqrt(ΓRe * (1. - fRe)) * ce,
    ]
    return L_ops
end
