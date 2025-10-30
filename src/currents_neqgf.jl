export current_particle_avg_neqgf
export current_heat_avg_neqgf

"""
Transmission function of the DQD according to eq. 13 of [Prech et. al. 2023]

# Arguments
- `ω::Real`: Energy (integration variable)

See also: [`build_dqd_leads`](@ref) and [`current_particle_avg_neqgf`](@ref)
"""
function _transmission_dqd(ω, dqd_leads::DqdLeads)
    ϵL, ϵR = get_onsite_energies(dqd_leads.dqd)

    zL = ω - ϵL + 0.5 * im * dqd_leads.ΓL
    zR = ω - ϵR + 0.5 * im * dqd_leads.ΓR

    numerator = dqd_leads.ΓL * dqd_leads.ΓR * dqd_leads.dqd.tc^2
    denominator = abs2(zL * zR - dqd_leads.dqd.tc^2)
    return numerator / denominator
end

"""
Particle current integrand for NEGFs approach
"""
function _integrand_particle_current_avg_neqgf(ω, dqd_leads::DqdLeads; left::Bool = true)
    μL, μR = get_chemical_potentials(dqd_leads.leads)

    Tω = _transmission_dqd(ω, dqd_leads)
    fL = fermi(ω, μL, dqd_leads.leads.TL)
    fR = fermi(ω, μR, dqd_leads.leads.TR)
    Δf = left ? fL - fR : fR - fL

    return Tω * Δf / (2π)  # sign depends on convention
end

@doc raw"""
Compute steady-state particle current through the non-interacting DQD using NEQGFs
according to eq. (14) [Prech et. al. 2023]

# Returns
- `I_avg::Real`
"""
function current_particle_avg_neqgf(dqd_leads::DqdLeads, left::Bool = true)
    integrand(ω) = _integrand_particle_current_avg_neqgf(ω, dqd_leads; left = left)
    I_avg, err = quadgk(integrand, -50, 50, atol=1e-8)
    if iszero(imag.(I_avg))
        return real.(I_avg)
    else
        error("Expected a real particle current, but got some complex value with nonzero imaginary part")
    end
end

"""
Particle current integrand for NEGFs approach
"""
function _integrand_heat_current_avg_neqgf(ω::Real, μ::Real, dqd_leads::DqdLeads; left::Bool = true)
    μL, μR = get_chemical_potentials(dqd_leads.leads)

    Tω = _transmission_dqd(ω, dqd_leads)
    fL = fermi(ω, μL, dqd_leads.leads.TL)
    fR = fermi(ω, μR, dqd_leads.leads.TR)
    Δf = left ? fL - fR : fR - fL

    return (ω - μ) * Tω * Δf / (2π)
end

@doc raw"""
Compute steady-state heat current through the non-interacting DQD using NEQGFs
according to eq. (B.3) [Potts et. al. 2021])

# Returns
- `I_avg::Real`
"""
function current_heat_avg_neqgf(dqd_leads::DqdLeads; left::Bool = true)
    μL, μR = get_chemical_potentials(dqd_leads.leads)
    μ = left ? μL : μR
    integrand_α(ω) = _integrand_heat_current_avg_neqgf(ω, μ, dqd_leads; left = left)

    J_avg, err = quadgk(integrand_α, -50, 50, atol=1e-8)

    if iszero(imag.(J_avg))
        return real.(J_avg)
    else
        error("Expected a real heat current, but got some complex value with nonzero imaginary part")
    end
end