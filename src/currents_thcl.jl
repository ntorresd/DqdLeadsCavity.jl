export current_particle_avg_thcl
export current_heat_avg_thcl

@doc raw"""
Calculates analytical steady-state particle current for the non-interacting local LME according to eq. (A12) [Prech et. al. 2023]

# Returns
- Steady-state average of the current operator

See also: [`build_dqd_basis`](@ref)
"""
function current_particle_avg_thcl(dqd_leads::DqdLeads; left::Bool = true)
    # Parameters L-R
    (; ϵL, ϵR, μL, μR, ΓL, ΓR, TL, TR) = context_LR(dqd_leads)
    ϵ_avg = dqd_leads.dqd.ϵ_avg
    fL = fermi(ϵ_avg, μL, TL)
    fR = fermi(ϵ_avg, μR, TR)
    Δf = left ? fL - fR : fR - fL

    δ = 2. * dqd_leads.dqd.Δϵ / (ΓL + ΓR) 
    numerator = 4 * dqd_leads.dqd.tc^2 * (ΓL * ΓR)^2 * Δf
    denominator = (ΓL + ΓR) * (4 * dqd_leads.dqd.tc^2 + (ΓL * ΓR) * (δ^2 + 1.))
    return numerator / denominator
end

@doc raw"""
Calculates analytical steady-state heat current through the non-interacting DQD
for the thermodynamically consistent local LME according to eq. (B.6) [Potts et. al. 2021]

# Returns
- Steady-state average of the current operator

See also: [`build_dqd_basis`](@ref)
"""
function current_heat_avg_thcl(dqd_leads::DqdLeads; left::Bool = true)
    # Parameters L-R
    (; ϵL, ϵR, μL, μR, ΓL, ΓR, TL, TR) = context_LR(dqd_leads)
    ϵ_avg = dqd_leads.dqd.ϵ_avg
    tc = dqd_leads.dqd.tc

    fL = fermi(ϵ_avg, μL, TL)
    fR = fermi(ϵ_avg, μR, TR)
    Δf = left ? fL - fR : fR - fL
    ϵα = left ? ϵ_avg - μL : ϵ_avg - μR

    Γ_avg = (ΓL + ΓR) / 2.
    γ = sqrt(ΓL * ΓR / Γ_avg^2) # Adimensional constant

    numerator = 2 * Γ_avg * (γ * tc)^2
    denominator = 4 * tc^2 + (γ * Γ_avg)^2 + (γ * dqd_leads.dqd.Δϵ)^2

    Jα = ϵα * Δf * numerator / denominator
    return Jα
end