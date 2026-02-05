using CairoMakie, LaTeXStrings

begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"\dot{N}",
        ylabel = L"I_L",
        title = "Particle current - Coherent drive",
        xscale = log10
    )
    lines!(ax, Ndot_range, I_L_num; color=:blue, linestyle=:solid, alpha=0.5, label=L"I_L^{num}")
    lines!(ax, Ndot_range, I_L_ld_num; color=:blue, linestyle=:dashdot, alpha=0.7, label=L"I_L^{num, ld}")
    lines!(ax, Ndot_range, I_L_ld_ana; color=:black, linestyle=:dash, alpha=0.7, label=L"I_L^{ana, ld}")
    axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "I_L_driven.png"), f)
    else
        display(f)
    end
end
