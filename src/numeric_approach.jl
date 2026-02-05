export get_particle_current_num
export get_heat_current_num
export current_output # not tested

@doc raw"""
Numerical particle current calculation
"""
function get_particle_current_num(
    dqd_leads::DqdLeads,
    ρ::QuantumObject{Operator},
    L_ops::AbstractVector{<:QuantumObject{Operator}}
)
    # operators
    dL, dR = build_dqd_fermi_ops_LR(dqd_leads.dqd)
    N_op = dL' * dL + dR' * dR
    # compute the particle current
    I = 0.0
    for L_op in L_ops
        I += expect(N_op, D_sop(L_op, ρ))
    end
    # check whether the current is real
    return _is_real(I)
end
function get_particle_current_num(
    dqd_leads_cavity::DqdLeadsCavityObj,
    ρ::QuantumObject{Operator},
    L_ops::AbstractVector{<:QuantumObject{Operator}}
)
    # operators
    dL, dR = build_dqd_fermi_ops_LR(dqd_leads_cavity)
    N_op = dL' * dL + dR' * dR
    # compute the particle current
    I = 0.0
    for L_op in L_ops
        I += expect(N_op, D_sop(L_op, ρ))
    end
    # check whether the current is real
    return _is_real(I)
end

@doc raw"""
Numerical heat current calculation
"""
function get_heat_current_num(
    dqd_leads::DqdLeads,
    ρ::QuantumObject{Operator},
    L_ops::AbstractVector{<:QuantumObject{Operator}},
    H_td::QuantumObject{Operator};
    left::Bool = true
)
    # parameters
    μL, μR = get_chemical_potentials(dqd_leads.leads)
    μ = left ? μL : μR
    # operators
    dL, dR = build_dqd_fermi_ops_LR(dqd_leads.dqd)
    N_op = dL' * dL + dR' * dR
    # compute the heat current
    J = 0.0
    for L_op in L_ops
        J += expect(H_td - μ * N_op, D_sop(L_op, ρ))
    end
    # check whether the current is real
    return _is_real(J)
end
function get_heat_current_num(
    dqd_leads_cavity::DqdLeadsCavityObj,
    ρ::QuantumObject{Operator},
    L_ops::AbstractVector{<:QuantumObject{Operator}},
    H_td::QuantumObject{Operator};
    left::Bool = true
)
    # parameters
    μL, μR = get_chemical_potentials(dqd_leads_cavity.dqd_leads.leads)
    μ = left ? μL : μR
    # operators
    dL, dR = build_dqd_fermi_ops_LR(dqd_leads_cavity)
    N_op = dL' * dL + dR' * dR
    # compute the heat current
    J = 0.0
    for L_op in L_ops
        J += expect(H_td - μ * N_op, D_sop(L_op, ρ))
    end
    # check whether the current is real
    return _is_real(J)
end

export current_output

@doc raw"""
Calculate the average of the output stochastic current according to eq. (40) of [Landi et. al. 2024]

# Arguments
- `ρ::QuantumObject{Operator}`: Reduced density matrix of the system
- `L_ops::Vector{QuantumObject{Operator}}}`: Vector of Lindblad operators

# Keyword Arguments
- `ν_vec::Vector{Int64}`: Weights corresponding to each Lindblad operator in `L_ops`

See also: [`build_L_ops_local`](@ref), [`build_L_ops_semilocal`](@ref)

# Returns
- `I_avg::Float64`: Output average stochastic current 
"""
function current_output(
    ρ::QuantumObject{Operator},
    L_ops::Vector{QuantumObject};
    ν_vec::Vector{Int}
)
    I_avg = 0im
    for i in eachindex(ν_vec)
        I_avg = I_avg + ν_vec[i] * expect(L_ops[i]' * L_ops[i], ρ)
    end
    if iszero(imag.(I_avg))
        return real.(I_avg)
    else
        error("Expected a real current, but got some complex value with nonzero imaginary part: $(I_avg)")
    end
end
