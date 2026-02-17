using CairoMakie, LaTeXStrings

begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"Γ",
        ylabel = L"\left< d_l^\dagger d_l' \right>",
        title = L"\text{DQD populations (local) vs } Γ = Γ_L / Γ_R",
        yscale = log10,
        xscale = log10
    )
    lines!(ax, Γ_range, n0_num_loc; label=L"\bar{n}_0^{num}", color=:green, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, n0_ana_loc; label=L"\bar{n}_0^{ana}", color=:green, linestyle=:dash, alpha=0.7)
	lines!(ax, Γ_range, nL_num_loc; label=L"\bar{n}_L^{num}", color=:blue, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, nL_ana_loc; label=L"\bar{n}_L^{ana}", color=:blue, linestyle=:dash, alpha=0.7)
	lines!(ax, Γ_range, nR_num_loc; label=L"\bar{n}_R^{num}", color=:red, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, nR_ana_loc; label=L"\bar{n}_R^{ana}", color=:red, linestyle=:dash, alpha=0.7)
	lines!(ax, Γ_range, αLR_num_loc; label=L"|α_{LR}^{num}|", color=:purple, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, αLR_ana_loc; label=L"|α_{LR}^{ana}|", color=:purple, linestyle=:dash, alpha=0.7)
    lines!(ax, Γ_range, n0_num_loc + nL_num_loc + nR_num_loc; color=:black, linestyle=:dash)
	axislegend(ax, position=:lb)
    if save_fig
        save(joinpath(@__DIR__, "pop_LR_loc.png"), f)
    else
        display(f)
    end
end

begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"Γ",
        ylabel = L"\left< d_l^\dagger d_l' \right>",
        title = L"\text{DQD populations (odd) vs } Γ = Γ_L / Γ_R",
        # yscale = log10,
        xscale = log10
    )
	lines!(ax, Γ_range, n0_num_odd; label=L"\bar{n}_0^{odd}", color=:green, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, nL_num_odd; label=L"\bar{n}_L^{odd}", color=:blue, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, nR_num_odd; label=L"\bar{n}_R^{odd}", color=:red, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, αLR_num_odd; label=L"|α_{LR}^{odd}|", color=:purple, linestyle=:solid, alpha=0.7)
    lines!(ax, Γ_range, n0_num_odd + nL_num_odd + nR_num_odd; color=:black, linestyle=:dash)
	axislegend(ax, position=:rt)
    if save_fig
        save(joinpath(@__DIR__, "pop_LR_odd.png"), f)
    else
        display(f)
    end
end

begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"Γ",
        ylabel = L"\left< d_l^\dagger d_l' \right>",
        title = L"\text{DQD populations vs } Γ = Γ_L / Γ_R",
        xscale = log10
    )
	lines!(ax, Γ_range, nL_num_odd; label=L"\bar{n}_L^{odd}", color=:blue, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, nL_ana_loc; label=L"\bar{n}_L^{ana}", color=:blue, linestyle=:dash, alpha=0.7)
	lines!(ax, Γ_range, nR_num_odd; label=L"\bar{n}_R^{odd}", color=:red, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, nR_ana_loc; label=L"\bar{n}_R^{ana}", color=:red, linestyle=:dash, alpha=0.7)
	lines!(ax, Γ_range, αLR_num_odd; label=L"|α_{LR}^{odd}|", color=:purple, linestyle=:solid, alpha=0.7)
	lines!(ax, Γ_range, αLR_ana_loc; label=L"|α_{LR}^{ana}|", color=:purple, linestyle=:dash, alpha=0.7)
	axislegend(ax, position=:rt)
    if save_fig
        save(joinpath(@__DIR__, "pop_LR_sloc_odd.png"), f)
    else
        display(f)
    end
end
