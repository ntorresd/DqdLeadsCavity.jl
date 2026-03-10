export get_particle_current_num
export get_heat_current_num
export current_output # not tested
export get_populations_LR_num, get_coherence_LR_num, get_concurrence_LR_num

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
        I += expect(N_op, D_sop(L_op) * ρ)
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
        I += expect(N_op, D_sop(L_op) * ρ)
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
        J += expect(H_td - μ * N_op, D_sop(L_op) * ρ)
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
        J += expect(H_td - μ * N_op, D_sop(L_op) * ρ)
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


function get_fermi_ge(dqd_leads::DqdLeads)
    ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
    μL, μR = get_chemical_potentials(dqd_leads)
    TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR

    # fLg, fLe, fRg, fRe
    return fermi(ϵg, μL, TL), fermi(ϵe, μL, TL), fermi(ϵg, μR, TR), fermi(ϵe, μR, TR)
end

@doc raw"""
Numerical populations of the DQD in the left-right basis
"""
function get_populations_LR_num(
	dqd_leads::DqdLeads,
	ρ::QuantumObject{Operator}
)
	ket_0, ket_L, ket_R, ket_D = build_dqd_basis_LR(dqd_leads.dqd)
	
	p0 = abs(ket_0' * ρ * ket_0)
	pL = abs(ket_L' * ρ * ket_L)
	pR = abs(ket_R' * ρ * ket_R)
	pd = abs(ket_D' * ρ * ket_D)

	return p0, pL, pR, pd
end

@doc raw"""
Numerical coherence between left-right states of the DQD
"""
function get_coherence_LR_num(
	dqd_leads::DqdLeads,
	ρ::QuantumObject{Operator}
)
	ket_0, ket_L, ket_R, ket_D = build_dqd_basis_LR(dqd_leads.dqd)
	α = abs(ket_L' * ρ * ket_R)

	return α

end

@doc raw"""
Numerical concurrence in the 
"""
function get_concurrence_LR_num(
	dqd_leads::DqdLeads,
	ρ::QuantumObject{Operator}
)
	p0, pL, pR, pd = get_populations_LR_num(dqd_leads, ρ)
	α = get_coherence_LR_num(dqd_leads, ρ)
	C = 2 * max(0, α - sqrt(p0 * pd))

	return C
end