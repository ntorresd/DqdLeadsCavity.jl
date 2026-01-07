using QuantumToolbox, DQDTransport
using Plots, LaTeXStrings


# include("../params/zenelaj2022.jl")
include("../params/prech2023.jl")
# --- local operators ---

# basis {|0>, |L>, {R>}
empty, left, right = basis(3, 0), basis(3, 1), basis(3, 2);

# population operators
nL = left*left'
nR = right*right'
# the only allowed transitions are 0 <-> L and 0 <-> R
sL = empty*left'   # |0><L|
sR = empty*right'  # |0><R|


"H and c_ops for DQD + cavity (displaced field)"
function H_cops_dqd(
    TL::Real, ΓL::Real,
    TR::Real, ΓR::Real,
    ϵ::Real, ϵ_avg::Real, Δμ::Real,
    tc::Real
)     
    Ω = sqrt(4*tc^2 + ϵ^2);
    θ = acos(-ϵ/Ω);

    sg = cos(θ/2)*sL - sin(θ/2)*sR;   # |0><g|
    se = sin(θ/2)*sL + cos(θ/2)*sR;   # |0><e|
    σz = se'*se - sg'*sg;

    # Hamiltonian
    Htot = Ω/2 * σz;

    # DQD
    Eg = -Ω/2 + ϵ_avg;
    Ee = Ω/2 + ϵ_avg;

    fLg = fermi(Eg, Δμ/2., TL);
    fRg = fermi(Eg, - Δμ/2., TR);

    fLe = fermi(Ee, Δμ/2, TL);
    fRe = fermi(Ee, - Δμ/2., TR);
    
    Γg0 = ΓL * fLg * cos(θ/2)^2 + ΓR * fRg * sin(θ/2)^2;
    Γe0 = ΓL * fLe * sin(θ/2)^2 + ΓR * fRe * cos(θ/2)^2;
    Γ0g = ΓL * (1. - fLg) * cos(θ/2)^2 + ΓR * (1. - fRg) * sin(θ/2)^2; 
    Γ0e = ΓL * (1. - fLe) * sin(θ/2)^2 + ΓR * (1. - fRe) * cos(θ/2)^2; 

    c_ops = [sqrt(Γg0)*sg', sqrt(Γ0g)*sg, sqrt(Γ0e)*se, sqrt(Γe0)*se'];
    return Htot, c_ops
end


"current from the left contact for DQD + cavity"
function current_L_dqd(
    ρ::QuantumObject,
    TL::Real, ΓL::Real,
    TR::Real, ΓR::Real,
    ϵ::Real, ϵ_avg::Real, Δμ::Real,
    tc::Real
)

    Ω = sqrt(4*tc^2 + ϵ^2);
    θ = acos(-ϵ/Ω);

    sg = cos(θ/2)*sL - sin(θ/2)*sR;   # |0><g|
    se = sin(θ/2)*sL + cos(θ/2)*sR;   # |0><e|

    Eg = -Ω/2 + ϵ_avg;
    Ee = Ω/2 + ϵ_avg;

    fLg = fermi(Eg, Δμ/2., TL);
    fRg = fermi(Eg, - Δμ/2., TR);

    fLe = fermi(Ee, μL, TL);
    fRe = fermi(Ee, μR, TR);

    Γg0L = ΓL * fLg * cos(θ/2)^2
    Γe0L = ΓL * fLe * sin(θ/2)^2
    Γ0gL = ΓL * (1. - fLg) * cos(θ/2)^2
    Γ0eL = ΓL * (1. - fLe) * sin(θ/2)^2
    
    Lrho = Γg0L * dissipator(sg', ρ) + Γ0gL * dissipator(sg, ρ) + Γe0L * dissipator(se', ρ) + Γ0eL * dissipator(se, ρ);
    real(expect(nL + nR, Lrho))
end

N_range = 10000
tc_range = range(1e-5, 10, N_range)

ρ_ss_list = [steadystate(
    H_cops_dqd(
        TL, ΓL,
        TR, ΓR,
        ϵ, ϵ_avg, Δμ,
        tc
    )...) for tc in tc_range
];

current_L = [current_L_dqd(
        ρ_ss_list[i],
        TL, ΓL,
        TR, ΓR,
        ϵ, ϵ_avg, Δμ,
        tc_range[i]
    ) for i in eachindex(ρ_ss_list)
]

plot(tc_range, current_L, xaxis = :log)

#--- Stability diagrams ---
N_squares = 250

x_max = 10 * T
y_max = 10 * T
x_range = range(- x_max, x_max, N_squares)
y_range = range(- y_max, y_max, N_squares)

" Get populations"
function get_populations(
    TL, ΓL,
    TR, ΓR,
    ϵ, ϵ_avg, Δμ,
    tc
)
    H_dqd, L_ops = H_cops_dqd(
        TL, ΓL,
        TR, ΓR,
        ϵ, ϵ_avg, Δμ,
        tc
    )
    return get_populations_LR(
        steadystate(H_dqd, L_ops);
        blockade = true
    )
end


populations_matrix = [
    get_populations(
        TL, ΓL,
        TR, ΓR,
        ϵL - ϵR, (ϵL + ϵR) / 2., Δμ,
        tc
    ) for ϵL in x_range, ϵR in y_range
]

heatmaps_plot = plot_populations_heatmaps(
    populations_matrix;
    x_range = x_range, y_range = y_range,
    xlabel = L"ϵ_L", ylabel = L"ϵ_R"
);
plot!(
    heatmaps_plot,
    size = (800, 700)
)

savefig(
    heatmaps_plot,
    "./external/plots/sd_ge_basis_prech2023.png"
)

# --- Current ---