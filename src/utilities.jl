export fermi
export id_no_vacuum
export context_LR
export D_sop

@doc raw"""
Fermi distribution

# Arguments
- `ϵ::Real`: Energy
- `μ::Real`: Chemical potential
- `T::Real`: Temperature
"""
function fermi(ϵ, μ, T)
    return 1 / (1 + exp((ϵ - μ) / T))
end

@doc raw"""
Lindblad dissipator superoperator
"""
function D_sop(x::QuantumObject{Operator}, ρ::QuantumObject{Operator})
    return x * ρ * x' - 1 / 2 * x' * x * ρ - 1 / 2 * ρ * x' * x
end

@doc raw"""
Identity for the subspace without the vacuum state
"""
function id_no_vacuum(dim::Int)
    return QuantumObject(Matrix(Diagonal(vcat(0.0, ones(dim - 1)))))
end

@doc raw"""
Helper function to gather g-e labelled quantities
"""
@inline function context_ge(dqd_leads::DqdLeads)
    ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
    ΓLg, ΓLe, ΓRg, ΓRe = get_transition_rates_ge(dqd_leads)
    fLg, fLe, fRg, fRe = get_fermi_ge(dqd_leads)
    (; ϵg, ϵe, ΓLg, ΓLe, ΓRg, ΓRe, fLg, fLe, fRg, fRe)
end

"""
Helper function to check whether a number is approximately real
"""
function _is_real(x::ComplexF64; atol::Float64 = 1e-20)
    if isapprox(imag(x), 0.0; atol = atol)
        return real(x)
    else
        error("Expected a real value, but got some complex value with nonzero imaginary part: $(x)")
    end
end
