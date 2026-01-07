using Revise
using DqdLeadsCavity, QuantumToolbox 
using Plots, LaTeXStrings

# --- Parameters SD ---
Γ = 1.
# DQD
ϵ_avg = 0.
tc = 1. * Γ
# U = 1e3 * tc
# U = 50. * tc
U = 100 * Γ
# Leads
# Δμ = 50. * tc
Δμ = 0.
# DQD-Leads
# Γ = 1.
# T = 10. * tc
T = 1. * Γ

# ---- L-R experimental ----
"""
Experimental Lindblad operators for the interacting DQD in the L-R basis
"""
function build_L_ops_nthcl_int(dqd_leads::DqdLeads)
    # Parameters g-e
    (; ϵL, ϵR, μL, μR, ΓL, ΓR, TL, TR) = context_LR(dqd_leads)
    μL, μR = get_chemical_potentials(dqd_leads)
    TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR
    # Fermi Functions
    fL = fermi(ϵL, μL, TL)
    fR = fermi(ϵR, μR, TR)
    fL_U = fermi(ϵL + U / 2., μL, TL)
    fR_U = fermi(ϵR + U / 2., μR, TR)    

    # Creation/Annihilation ops (for each dot)
    cL, cR = build_dqd_ladder_ops_LR(dqd_leads.dqd)

    # Lindblad jump operators (lead-dot tunneling)
    L_ops = [
        # Left lead jump operators
        sqrt(ΓL * fL) * (1. - cR' * cR) * cL',
        sqrt(ΓL * (1. - fL)) * (1. - cR' * cR) * cL,
        sqrt(ΓL * fL_U) * cR' * cR * cL',
        sqrt(ΓL * (1. - fL_U)) * cR' * cR * cL,
        ## Right Lead jump operators
        sqrt(ΓR * fR) * (1. - cL' * cL) * cR',
        sqrt(ΓR * (1. - fR)) * (1. - cL' * cL) * cR,
        sqrt(ΓR * fR_U) * cL' * cL * cR',
        sqrt(ΓR * (1. - fR_U)) * cL' * cL * cR,
    ]

    return L_ops
end

function get_populations_LR(
    ρ::QuantumObject{Operator},
    dqd::Dqd
)
    nL, nR, nD = build_dqd_number_ops_LR(dqd)

    populations = real.([
        1 - expect(nL + nR - nD, ρ),    # p0
        expect(nL - nD, ρ),             # pL
        expect(nR - nD, ρ),             # pR
        expect(nD, ρ)                   # pD
    ])

    return Tuple(populations)
end

# ---- Stability Diagrams ----
N_squares = 250
max_range = 100 * tc
# max_range = 100 * T
x_range = range(-max_range, max_range, N_squares)
y_range = range(-max_range, max_range, N_squares)
# x_range = [1.]
# y_range = [-0.3]

function get_populations_LR(leads::Leads; ϵL::Real, ϵR::Real)
    dqd = Dqd(ϵL - ϵR, (ϵL + ϵR) / 2., tc, U, false)
    dqd_leads = DqdLeads(dqd, leads, Γ, Γ)
    # show(dqd_leads)
    H_dqd = build_H_dqd_LR(dqd)
    L_ops = build_L_ops_nthcl_int(dqd_leads)
    ρ = steadystate(H_dqd, L_ops)

    nL, nR, nD = build_dqd_number_ops_LR(dqd_leads.dqd)
    p0 = real(1 - expect(nL + nR - nD, ρ))
    pL = real(expect(nL - nD, ρ))
    pR = real(expect(nR - nD, ρ))
    pD = real(expect(nD, ρ))

    return p0, pL, pR, pD
end

leads = Leads(T, T, Δμ)
populations_matrix = [get_populations_LR(leads; ϵL = ϵL, ϵR = ϵR) for ϵL in x_range, ϵR in y_range]

# ---- Plot populations heatmaps ----
# Extract populations
p0 = getindex.(populations_matrix, 1);
pL = getindex.(populations_matrix, 2);
pR = getindex.(populations_matrix, 3);
pD = getindex.(populations_matrix, 4);
p_max = maximum(maximum.([p0, pL, pR, pD]));
# Plot populations
heat_maps = plot(
    heatmap(x_range, y_range, p0', title = L"|0\rangle"),
    heatmap(x_range, y_range, pD', title = L"|D\rangle"),
    heatmap(x_range, y_range, pL', title = L"|L\rangle"),
    heatmap(x_range, y_range, pR', title = L"|R\rangle"),
    layout = (2, 2),
    size = (800, 700),
    xlabel = L"\epsilon_L",
    ylabel = L"\epsilon_R",
    clim = (0, p_max)
)

