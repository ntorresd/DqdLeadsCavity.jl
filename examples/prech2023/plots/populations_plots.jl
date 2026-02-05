using CairoMakie, LaTeXStrings

begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"t_c",
        ylabel = L"\left< d_\sigma^\dagger d_\sigma' \right>",
        title = L"\text{DQD populations (global) vs } t_c",
        xscale = log10
    )
	lines!(ax, tc_range, n_dqd_num_g_gl; label=L"\bar{n}_g^{num}", color=:blue, linestyle=:solid, alpha=0.7)
	lines!(ax, tc_range, n_dqd_ana_g_gl; label=L"\bar{n}_g^{ana}", color=:blue, linestyle=:dash, alpha=0.7)
	lines!(ax, tc_range, n_dqd_num_e_gl; label=L"\bar{n}_e^{num}", color=:red, linestyle=:solid, alpha=0.7)
	lines!(ax, tc_range, n_dqd_ana_e_gl; label=L"\bar{n}_e^{ana}", color=:red, linestyle=:dash, alpha=0.7)
	lines!(ax, tc_range, α_ge_num_gl; label=L"|α_{ge}^{num}|", color=:purple, linestyle=:solid, alpha=0.7)
	lines!(ax, tc_range, α_ge_ana_gl; label=L"|α_{ge}^{ana}|", color=:purple, linestyle=:dash, alpha=0.7)
	axislegend(ax, position=:lt)
    if save_fig
        save(joinpath(@__DIR__, "pop_ge_gl.png"), f)
    else
        display(f)
    end
end

begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"t_c",
        ylabel = L"\left< d_l^\dagger d_l' \right>",
        title = L"\text{DQD populations (global) vs } t_c",
        xscale = log10
    )
	lines!(ax, tc_range, n_dqd_num_L_gl; label=L"\bar{n}_L^{num}", color=:blue, linestyle=:solid, alpha=0.7)
	lines!(ax, tc_range, n_dqd_ana_L_gl; label=L"\bar{n}_L^{ana}", color=:blue, linestyle=:dash, alpha=0.7)
	lines!(ax, tc_range, n_dqd_num_R_gl; label=L"\bar{n}_R^{num}", color=:red, linestyle=:solid, alpha=0.7)
	lines!(ax, tc_range, n_dqd_ana_R_gl; label=L"\bar{n}_R^{ana}", color=:red, linestyle=:dash, alpha=0.7)
	lines!(ax, tc_range, α_LR_num_gl; label=L"|α_{LR}^{num}|", color=:purple, linestyle=:solid, alpha=0.7)
	lines!(ax, tc_range, α_LR_ana_gl; label=L"|α_{LR}^{ana}|", color=:purple, linestyle=:dash, alpha=0.7)
	axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "pop_LR_gl.png"), f)
    else
        display(f)
    end
end
