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
Δϵ = 0.
ϵ_avg = 0.
tc = 6.5 * Γ

# --- Template system ----
dqd = Dqd(Δϵ, ϵ_avg, tc)
leads = Leads(T, T, Δμ)
dqd_leads = DqdLeads(dqd, leads, Γ, Γ)

# Heat current
function current_heat_tdbk(
    ρ::QuantumObject,
    H_bk::QuantumObject,
    N_op::QuantumObject,
    L_ops,
    μ_α::Real
)
    J_list = Vector{Float64}(undef, length(L_ops))
    for (i, L_op) in enumerate(L_ops)
        J_list[i] = real(expect((H_bk - μ_α * N_op) * L_op, ρ))
    end
    return J_list
end

# H_dqd = build_H_dqd_ge(dqd_leads.dqd)
# ng, ne, nD = build_dqd_number_ops_ge(dqd_leads.dqd)
# N_op = ng + ne
# L_ops = build_L_ops_thcg_int(dqd_leads)
# ρ = steadystate(H_dqd, L_ops)
# J_thcg = current_heat_tdbk(ρ, H_dqd, N_op, L_ops, μL)

N_range = Int(1e4)
J_thcg_L = Vector{Float64}(undef, N_range)
J_thcg_R = Vector{Float64}(undef, N_range)
x_range = logrange(1e-10, 1, N_range)
for (i, x) in enumerate(x_range)
    dqd_leads.dqd.tc = x_range[i]
    μL, μR = get_chemical_potentials(dqd_leads.leads)
    ng, ne, nD = build_dqd_number_ops_ge(dqd_leads.dqd)
    N_op = ng + ne
    H_dqd = build_H_dqd_ge(dqd_leads.dqd)
    L_thcg = build_L_ops_thcg_int(dqd_leads)
    ρss = steadystate(H_dqd, L_thcg)

    J_thcg = current_heat_tdbk(ρss, H_dqd, N_op, L_thcg, μL)
    J_thcg_L[i] = sum(J_thcg[1:8])
end


xlabel = L"t_c / \Gamma"
legend = :left
xaxis = :log
populations_plot = plot(
    xlabel = xlabel,
    ylabel = L"J",
    legend = legend,
    xaxis = xaxis,
    dpi = 200
);
plot!(
    x_range/Γ, J_thcg_L,
    label = L"J_L^{thcg}"
)






#---- Analytical Heat Currents ----
N_range = 10000
tc_range = logrange(1e-4, 25 * tc, N_range)
# Compute NEQF current (non-interacting)
J_neqfg = map(tc_range) do tc
    dqd.tc = tc
    current_heat_avg_neqgf(dqd_leads; left = false)
end
# Compute analytical heat current (local, non-interacting)
J_thcl = map(tc_range) do tc
    dqd.tc = tc
    current_heat_avg_thcl(dqd_leads; left = false)
end
# Compute analytical heat current (global, non-interacting)
J_thcg = map(tc_range) do tc
    dqd.tc = tc
    current_heat_avg_thcg_ana(dqd_leads; left = false)
end

current_heat_plot = plot(
	tc_range / Γ, [J_neqfg, J_thcl, J_thcg],
    xlabel = L"t_c / \Gamma",
    ylabel = L"\langle J \rangle / \Gamma",
	label = [L"J_{neqgf}" L"J_{thcl}" L"J_{thcg}"],
	linestyle = :dash,
    legend = :left,
    xaxis = :log,
    dpi = 200
)
plot!(
    [Δμ / 2.], seriestype = :vline,
    linecolor = :black, linestyle = :dash,
    label = ""
)


#---- Particle currents -----
N_range = 10000
tc_range = logrange(1e-4, 25 * tc, N_range)
# Compute NEQF current (non-interacting)
I_neqfg = map(tc_range) do tc
    dqd.tc = tc
    current_particle_avg_neqgf(dqd_leads)
end
# Compute analytical current (local, non-interacting)
I_thcl = map(tc_range) do tc
    dqd.tc = tc
    current_particle_avg_thcl(dqd_leads)
end

current_particle_plot = plot(
	tc_range / Γ, [I_neqfg, I_thcl],
    xlabel = L"t_c / \Gamma",
    ylabel = L"\langle I \rangle / \Gamma",
	label = [L"I_{neqgf}" L"I_{thcl}"],
	linestyle = :dash,
    legend = :topleft,
    xaxis = :log,
    dpi = 200
)
plot!(
    [Δμ / 2.], seriestype = :vline,
    linecolor = :black, linestyle = :dash,
    label = ""
)