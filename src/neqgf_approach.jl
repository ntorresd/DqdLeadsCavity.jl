export get_particle_current_neqgf
export get_current_heat_neqgf

"""
Transmission function of the DQD under the wideband approximation
[eq. (B5) Potts2021] ≡ [eq. (A49) Prech2023]

# Arguments
- `ω::Real`: Energy (integration variable)

See also: [`build_dqd_leads`](@ref) and [`get_particle_current_neqgf`](@ref)
"""
function _transmission_dqd(ω, dqd_leads::DqdLeads)
    ϵL, ϵR = get_onsite_energies(dqd_leads.dqd)
    ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
    tc = dqd_leads.dqd.tc

    zL = ω - ϵL + 0.5 * im * ΓL
    zR = ω - ϵR + 0.5 * im * ΓR

    numerator = ΓL * ΓR * tc^2
    denominator = abs2(zL * zR - tc^2)
    return numerator / denominator
end

"""
Particle current integrand for NEGFs approach
"""
function _integrand_particle_current_avg_neqgf(ω, dqd_leads::DqdLeads; left::Bool = true)
    μL, μR = get_chemical_potentials(dqd_leads)

    Tω = _transmission_dqd(ω, dqd_leads)
    fL = fermi(ω, μL, dqd_leads.leads.TL)
    fR = fermi(ω, μR, dqd_leads.leads.TR)
    Δf = left ? fL - fR : fR - fL

    return Tω * Δf / (2π)  # sign depends on convention
end

@doc raw"""
Compute steady-state particle current through the non-interacting DQD using NEQGFs
according to [eq. (14) Prech et. al. 2023]

# Returns
- `I_avg::Real`
"""
function get_particle_current_neqgf(
    dqd_leads::DqdLeads;
    int_lims::Tuple{<:Real,<:Real},
    atol::Float64 = 1e-8)
    # integrand functions
    integrand_L(ω) = _integrand_particle_current_avg_neqgf(ω, dqd_leads; left = true)
    integrand_R(ω) = _integrand_particle_current_avg_neqgf(ω, dqd_leads; left = false)

    ωmin, ωmax = int_lims
    if ((ωmax > 0.) && (ωmin == -ωmax)) == false
        error("Integration limits must be symmetric around zero, i.e., (−ωmax, ωmax)")
    end

    I_avg_L, err = quadgk(integrand_L, ωmin, ωmax; atol = atol)
    I_avg_R, err = quadgk(integrand_R, ωmin, ωmax; atol = atol)
    if iszero(imag.(I_avg_L)) && iszero(imag.(I_avg_R))
        return real.(I_avg_L), real.(I_avg_R)
    else
        error("Expected a real particle current, but got some complex value with nonzero imaginary part")
    end
end

"""
Particle current integrand for NEGFs approach
"""
function _integrand_heat_current_neqgf(ω::Real, μ::Real, dqd_leads::DqdLeads; left::Bool = true)
    μL, μR = get_chemical_potentials(dqd_leads)

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
function get_current_heat_neqgf(
    dqd_leads::DqdLeads;
    int_lims::Tuple{<:Real,<:Real},
    atol::Float64 = 1e-8
    )
    μL, μR = get_chemical_potentials(dqd_leads)

    integrand_L(ω) = _integrand_heat_current_neqgf(ω, μL, dqd_leads; left = true)
    integrand_R(ω) = _integrand_heat_current_neqgf(ω, μR, dqd_leads; left = false)

    ωmin, ωmax = int_lims
    if ((ωmax > 0.) && (ωmin == -ωmax)) == false
        error("Integration limits must be symmetric around zero, i.e., (−ωmax, ωmax)")
    end

    J_avg_L, err = quadgk(integrand_L, ωmin, ωmax; atol = atol)
    J_avg_R, err = quadgk(integrand_R, ωmin, ωmax; atol = atol)

    if iszero(imag.(J_avg_L)) && iszero(imag.(J_avg_R))
        return real.(J_avg_L), real.(J_avg_R)
    else
        error("Expected a real heat current, but got some complex value with nonzero imaginary part")
    end
end