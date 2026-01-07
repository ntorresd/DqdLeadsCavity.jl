using Revise
using DqdLeadsCavity, QuantumToolbox 
using Plots, LaTeXStrings

# --- Parameters [Prech et. al. 2023] ---
Γ = 1.
T = 10. * Γ
Δϵ = 0.         # DQD level splitting 
ϵ_avg = 0.      # Average energy wrt μ=0
Δμ = 6.5 * T    # Δμ ≡ eV
U = 1e6
tc = 1.

# --- Template system ----
dqd = Dqd(Δϵ, ϵ_avg, tc, U, false)
leads = Leads(T, T, Δμ)
dqd_leads = DqdLeads(dqd, leads, Γ, Γ)

# --- Populations ----

function get_populations_ge(dqd_leads::DqdLeads)
    H_dqd = build_H_dqd_ge(dqd_leads.dqd)
    L_ops = build_L_ops_thcg_int(dqd_leads)
    ρ = steadystate(H_dqd, L_ops)

    ng, ne, nD = build_dqd_number_ops_ge(dqd_leads.dqd)
    p0 = real(1 - expect(ng + ne - nD, ρ))
    pg = real(expect(ng - nD, ρ))
    pe = real(expect(ne - nD, ρ))
    pD = real(expect(nD, ρ))

    return p0, pg, pe, pD
end


# ---- Plot populations ----
N_range = 1000
populations_list = Vector{NTuple{4, Float64}}(undef, N_range)


PLOT_PATH = "pops-vs-tc_U=$U.png"
x_range = logrange(10, 3e2, N_range)
xlabel = L"t_c / \Gamma"
legend = :left
xaxis = :log
# PLOT_PATH = "pops-vs-U_tc=$tc.svg"
# x_range = range(0, 150, N_range)
# xlabel = L"U"
# legend = :left
# xaxis = :linear

populations_plot = plot(
    xlabel = xlabel,
    ylabel = L"p",
    legend = legend,
    xaxis = xaxis,
    dpi = 200
);

for (i, x) in enumerate(x_range)
    dqd_leads.dqd.tc = x
    # dqd_leads.dqd.U = x
    populations_list[i] = get_populations_ge(dqd_leads)
end


p0 = getindex.(populations_list, 1)
pg = getindex.(populations_list, 2)
pe = getindex.(populations_list, 3)
pd = getindex.(populations_list, 4)
p_max = maximum(maximum.([p0, pg, pe, pd]))

plot!(
    x_range, [p0, pg, pe, pd],
    label = [L"p_{0}" L"p_{g}" L"p_{e}" L"p_{d}"],
    linestyle = [:dash :solid :solid :dash],
    alpha = [0.7 0.7 0.7 0.7],
    dpi = 200
)

savefig(
    populations_plot,
    string("../plots/populations/thcg/", PLOT_PATH)
)
