using Revise, DqdLeadsCavity
using QuantumToolbox
using LinearAlgebra

run_plots = true
save_fig = false

# Setup
begin
	# DQD parameters
	local ϵ_avg = 0.0
	local tc = 1.
	# local Δϵ = - 17.57
    # local Δϵ = - 2 * tc
    local Δϵ = 0.0
    local blockade = true
	
	local γm = 0.			# Relaxation rate
	local γphi = 3.92 	# Dephasing rate
	
    # Leads
	local T = 0.01 * tc
	local Δμ = 6.5 * tc
	local μ_avg = 0.0
	local Γ = 0.1

	leads = Leads(T, T, Δμ, μ_avg)


	dqd = Dqd(Δϵ, ϵ_avg, tc, γm, γphi, blockade)

	# Cavity parameters (from [zenelaj2022])
	local ωc = get_Ω(dqd)
	local κ = 0.337
	local Tc = 1e-7
	# Coherent driving
	local κin = 0.094
	local Ndot = 0.01
	local ωd = get_Ω(dqd)
	local gJC = 0.457

	cavity = Cavity(ωc, κ, Tc, ωd, κin, Ndot, gJC, 3)
	
	dqd_leads = DqdLeads(dqd, leads, Γ, Γ)
	dqd_leads_cavity = DqdLeadsCavityObj(dqd_leads, cavity)
end

"""
Steady state analytical populations for the semi-local LME
according to [eq.(C5) Prech et. al. 2023].
This is only valid for ϵL=ϵR
"""
function populations_ana_semilocal(dqd_leads::DqdLeads)
    # parameters
    ϵ_avg = dqd_leads.dqd.ϵ_avg;
    tc = dqd_leads.dqd.tc;
	μL, μR = get_chemical_potentials(dqd_leads);
	TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR;
    U = dqd_leads.dqd.U;
	ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR;

    # Fermi distributions
    fL = fermi(ϵ_avg, μL, TL);
    fR = fermi(ϵ_avg, μR, TR);
    fL_U = fermi(ϵ_avg + U, μL, TL);
    fR_U = fermi(ϵ_avg + U, μR, TR);
    f = (ΓL * fL + ΓR * fR) / (ΓL + ΓR);
    f_U = (ΓL * fL_U + ΓR * fR_U) / (ΓL + ΓR);
    ΔL = fL - fL_U;
    ΔR = fR - fR_U;

    # auxiliar definitions
    a = ΓL * ΓR * (1 - f + f_U);
    b = (ΓL * fL_U * ΔR + ΓR * fR_U * ΔL) / (ΓL + ΓR);
    c = (fL_U * fR - fR_U * fL) / (ΓL + ΓR);
    d = 4. * tc^2 * f;
    g = ΓL * ΓR * (fL * (1. - fR_U) + fR * (1. - fL_U)) / (ΓL + ΓR);

    # populations 
    denominator = 4. * tc^2 * (1 + f - f_U) + ((ΓL * ΓR) * (1. - ΔL * ΔR) * (1. - f + f_U));

    pL = (d * (1. - f_U) + a * (fL * (1. - fR) + b + ΓR * c)) / denominator;
    pR = (d * (1. - f_U) + a * (fR * (1. - fL) - b + ΓL * c)) / denominator;
    pD = (d * f_U + a * ((ΓL * fL_U * fR + ΓR * fR_U * fL) / (ΓL + ΓR) - b)) / denominator;
    p0 = 1. - pL - pR - pD;
    α = 2. * im * tc * g / denominator;

    return p0, pL, pR, pD, α
end

"""
Unitary part of the Liouvillian
"""
function liouvillian_unit(
    H::QuantumObject{Operator};
    Id = I(prod(H.dimensions))
)
    return -im * (spre(H, Id,) - spost(H, Id))
end

"""
Lindblad dissipator superoperator for even sectors
"""
function D_sop_even(
    L_op::QuantumObject{Operator};
    Id = I(prod(L_op.dimensions))
)
    return 2. * sprepost(L_op, L_op') - (spre(L_op' * L_op, Id) + spost(L_op' * L_op, Id))
end

"""
Lindblad dissipator superoperator for odd sectors
"""
function D_sop_odd(
    L_op::QuantumObject{Operator};
    Id = I(prod(L_op.dimensions))
)
    return - 2. * sprepost(L_op, L_op') - (spre(L_op' * L_op, Id) + spost(L_op' * L_op, Id))
end

"""
Modified Liouvillian for odd contributions
"""
function liouvillian_odd(
    H::QuantumObject{Operator},
    L_ops::Union{Nothing,AbstractVector,Tuple};
    Id = I(prod(H.dimensions))
)
    L = liouvillian_unit(H; Id = Id)
    for L_op in L_ops
        L += D_sop_odd(L_op; Id = Id)
    end

    return L
end