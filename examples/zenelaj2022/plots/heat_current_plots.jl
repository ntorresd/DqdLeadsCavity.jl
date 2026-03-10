using CairoMakie, LaTeXStrings

begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"\dot{N}",
        ylabel = L"J",
        title = "Heat current vs Cavity drive",
        xscale = log10,
    )
    lines!(ax, Ndot_range, JL_num_drive; color=:blue, linestyle=:dash, alpha=0.7, label=L"J_L^{num}")
    lines!(ax, Ndot_range, JR_num_drive; color=:red, linestyle=:dash, alpha=0.7, label=L"J_R^{num}")
    lines!(ax, Ndot_range, J_ana_L_nd; color=:blue, linestyle=:solid, alpha=0.5, label=L"J_L^{ana, nd}")
    # lines!(ax, Ndot_range, J_ana_L_ld;    color=:blue, linestyle=:solid, alpha=0.5, label=L"J_L^{ana, ld}")
    axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "J_drive.png"), f)
    else
        display(f)
    end
end