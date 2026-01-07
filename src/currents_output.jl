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

# function current_heat_avg_thcg_num(
#     ρ::QuantumObject{Operator},
#     H_td::QuantumObject{Operator},
#     L_ops::Vector{QuantumObject{Operator}};
#     lindblad_indexes::Vector{Int}
# )
#     J = 0.0
#     for index in lindblad_indexes
#         J += real(expect(H_td * L_ops[index] * ρ))
#     end
# end