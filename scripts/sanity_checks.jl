using Revise
using DqdLeadsCavity, QuantumToolbox
using Plots, LaTeXStrings

# --- Parameters equal onsite energies ---
# Leads
Γ = 1.
T = 10 * Γ
TL = T
TR = T
Δμ =  2 * Γ

# Dqd
Δϵ = 0.5
ϵ_avg = 0.
tc = Γ

# --- Template system ----
dqd = Dqd(Δϵ, ϵ_avg, tc)
leads = Leads(T, T, Δμ)
dqd_leads = DqdLeads(dqd, leads, Γ, Γ)

# --- L/R basis ----
cL, cR =  build_dqd_ladder_ops_LR(dqd)
ket_0, ket_L, ket_R = build_dqd_basis_LR(dqd)
ket_v = cL * ket_0

(cL' * ket_0 == ket_L)
(cR' * ket_0 == ket_R)
(cL * ket_L == ket_0)
(cR * ket_R == ket_0)
(cR * ket_L == ket_v)
(cL * ket_R == ket_v)

# ---- g/e basis ----
cg, ce = build_dqd_ladder_ops_ge(dqd)
ket_0, ket_g, ket_e = build_dqd_basis_ge(dqd)
ket_v = cg * ket_g

ce * ket_e == ket_g

