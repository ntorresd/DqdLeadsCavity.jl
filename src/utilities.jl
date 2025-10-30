export fermi
export get_Δd, get_dim
export id_no_vacuum
export context_LR

@doc raw"""
Fermi distribution

# Arguments
- `ϵ::Real`: Energy
- `μ::Real`: Chemical potential
- `T::Real`: Temperature
"""
function fermi(ϵ, μ, T)
    return 1/(1 + exp((ϵ - μ) / T))
end

@doc raw"""
Detuning between the DQD and cavity drive
"""
function get_Δd(dqd::Dqd, cavity::Cavity)
    return get_Ω(dqd) - cavity.ωd
end

@doc raw"""
Identity for the subspace without the vacuum state
"""
function id_no_vacuum(dim::Int)
    return QuantumObject(Matrix(Diagonal(vcat(0.0, ones(dim - 1))))) 
end

@doc raw"""
Get object dimension
"""
function get_dim(dqd::Dqd)
    dim = dqd.blockade ? 3 : 4
    return dim
end

@doc raw"""
Helper function to gather L-R labelled quantities
"""
@inline function context_LR(dqd_leads::DqdLeads)
    ϵL, ϵR = get_onsite_energies(dqd_leads.dqd)
    μL, μR = get_chemical_potentials(dqd_leads.leads)
    (; ϵL, ϵR, μL, μR,
       ΓL = dqd_leads.ΓL, ΓR = dqd_leads.ΓR,
       TL = dqd_leads.leads.TL, TR = dqd_leads.leads.TR
    )
end
