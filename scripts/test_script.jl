using Revise
using DqdLeadsCavity, QuantumToolbox
using Plots, LaTeXStrings

# Prech2023

# Leads
Γ = 1.
T = 10. * Γ
TL = T
TR = T
Δμ = 6.5 * T

# Dqd
Δϵ = 0.1 * Γ
ϵ_avg = 0.
tc = 6.5

dqd = Dqd(Δϵ, ϵ_avg, tc, false);
leads = Leads(TL, TR, Δμ);
dqd_leads = DqdLeads(dqd, leads, Γ, Γ)

μL, μR = get_chemical_potentials(dqd_leads)

N_range = 10000;
tc_range = logrange(1e-4, 25 * tc, N_range);


function current_heat_tdbk(
    dqd_leads::DqdLeads, # TODO: (refac) Remove this from build_dqd_basis (dqd -> blockade)
    ρ::QuantumObject,
    H_tdbk::QuantumObject,
    L_op_α::QuantumObject,
    μ_α
)
    N_dqd = build_dqd_number_op(dqd_leads.dqd)
    return real(expect((H_tdbk - μ_α * N_dqd) * L_op_α, ρ))
end

J_thcg_num = map(tc_range) do tc
    dqd_leads.dqd.tc = tc
    H_dqd = build_H_dqd_ge(dqd_leads.dqd)
    L_thcg = build_L_ops_thcg(dqd_leads)
    ρss = steadystate(H_dqd, L_thcg)
    # print(ρss, "\n")
    current_heat_tdbk(dqd_leads, ρss, H_dqd, L_thcg[1] + L_thcg[2], μL)
end

plot(
    tc_range / Γ, J_thcg_num,
    xlabel = L"t_c / \Gamma",
    ylabel = L"\langle I \rangle / \Gamma",
    label = L"J_{thcg}^{num}",
    legend = :left,
    xaxis = :log,
    dpi = 200
)
