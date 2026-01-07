using Revise
using DqdLeadsCavity, QuantumToolbox
using Plots, LaTeXStrings

Δϵ = -17.57
ϵ_avg = 0.
tc = 6.78
Γ = 1.

Δμ = 0.
μ_avg = 0.
T = 1e-7

dqd = Dqd(Δϵ, ϵ_avg, tc)
leads = Leads(T, T, Δμ, μ_avg)
dqd_leads = DqdLeads(dqd, leads, Γ, Γ)

μL, μR = get_chemical_potentials(dqd_leads)

# --- Local operators ---
ket_0, ket_L, ket_R = build_dqd_basis_LR(dqd)
cL, cR = build_dqd_ladder_ops_LR(dqd)
nL, nR = build_dqd_number_ops_LR(dqd)

Ω = get_Ω(dqd)
θ = get_θ(dqd)

cg, ce = build_dqd_ladder_ops_ge(dqd)

sg = cos(θ/2)*sL0 - sin(θ/2)*sR0;   # |0><g|
se = sin(θ/2)*sL0 + cos(θ/2)*sR0;   # |0><e|