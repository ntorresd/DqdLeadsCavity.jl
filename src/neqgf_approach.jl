export get_particle_current_neqgf
export get_heat_current_neqgf
export get_coherence_LR_neqgf

"""
Transmission function of the DQD under the wideband approximation
[eq. (B5) Potts2021] ≡ [eq. (A49) Prech2023]

# Arguments
- `ω::Real`: Energy (integration variable)
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
Steady-state particle current through the non-interacting DQD using NEQGFs
according to [eq. (14) Prech et. al. 2023]

# Returns
- `I_avg::Real`
"""
function get_particle_current_neqgf(
    dqd_leads::DqdLeads;
    left::Bool = true,
    int_lims::Tuple{<:Real,<:Real},
    atol::Float64 = 1e-8
)
    # integrand functions
    integrand(ω) = _integrand_particle_current_avg_neqgf(ω, dqd_leads; left = left)

    ωmin, ωmax = int_lims
    if ((ωmax > 0.) && (ωmin == -ωmax)) == false
        error("Integration limits must be symmetric around zero, i.e., (−ωmax, ωmax)")
    end
    I, err = quadgk(integrand, ωmin, ωmax; atol = atol)
    return I
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
Steady-state heat current through the non-interacting DQD using NEQGFs
according to eq. (B.3) [Potts et. al. 2021])

# Returns
- `I_avg::Real`
"""
function get_heat_current_neqgf(
    dqd_leads::DqdLeads;
    left::Bool = true,
    atol::Float64 = 1e-8,
    int_lims::Tuple{<:Real,<:Real}
    )
    μL, μR = get_chemical_potentials(dqd_leads)
    μ = left ? μL : μR
    integrand(ω) = _integrand_heat_current_neqgf(ω, μ, dqd_leads; left = left)

    ωmin, ωmax = int_lims
    if ((ωmax > 0.) && (ωmin == -ωmax)) == false
        error("Integration limits must be symmetric around zero, i.e., (−ωmax, ωmax)")
    end

    J, err = quadgk(integrand, ωmin, ωmax; atol = atol)
    return J
end

function _integrand_coherence_LR_neqgf(ω, dqd_leads::DqdLeads)
    # parameters
    tc = dqd_leads.dqd.tc
    ϵL, ϵR = get_onsite_energies(dqd_leads.dqd)
    μL, μR = get_chemical_potentials(dqd_leads)
    TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR
    ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
    # derived parameters
    fL, fR = fermi(ϵL, μL, TL), fermi(ϵR, μR, TR)
    # lesser Green's function
    numerator = im * tc * ((ω - ϵR + im * ΓR / 2.) * ΓL * fL + (ω - ϵL - im * ΓL / 2.) * ΓR * fR)
    denominator = abs((ω - ϵL + im * ΓL / 2.) * (ω - ϵR + im * ΓR / 2.) - tc^2)^2
    G_lesser_LR = numerator / denominator

    return G_lesser_LR / (2. * pi)
end

@doc raw"""
Steady-state coherence between the left and right dot using NEQGFs
according to eq. (A.43) [Prech2023])
"""
function get_coherence_LR_neqgf(
    dqd_leads::DqdLeads;
    atol::Float64 = 1e-8,
    int_lims::Tuple{<:Real,<:Real}
)
    integrand(ω) = _integrand_coherence_LR_neqgf(ω, dqd_leads)

    ωmin, ωmax = int_lims
    if ((ωmax > 0.) && (ωmin == -ωmax)) == false
        error("Integration limits must be symmetric around zero, i.e., (−ωmax, ωmax)")
    end

    α, err = quadgk(integrand, ωmin, ωmax; atol = atol)
    return abs(α)
end