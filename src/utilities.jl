export fermi
export get_Δd, get_dim
export id_no_vacuum

"""
Fermi distribution

# Arguments
- `ϵ::Real`: Energy
- `μ::Real`: Chemical potential
- `T::Real`: Temperature
"""
function fermi(ϵ, μ, T)
    return 1/(1 + exp((ϵ - μ) / T))
end

"""
Detuning between the DQD and cavity drive
"""
function get_Δd(dqd::Dqd, cavity::Cavity)
    return get_Ω(dqd) - cavity.ωd
end

"""
Identity for the subspace without the vacuum state
"""
function id_no_vacuum(dim::Int)
    return QuantumObject(Matrix(Diagonal(vcat(0.0, ones(dim - 1))))) 
end

"""
Get object dimension
"""
function get_dim(dqd::Dqd)
    dim = dqd.blockade ? 3 : 4
    return dim
end