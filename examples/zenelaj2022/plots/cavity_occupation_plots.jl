using CairoMakie, LaTeXStrings

begin
	f = Figure(size = (900, 600))
	ax = Axis(
		f[1, 1];
		xlabel = L"\dot{N}",
		ylabel = L"\left< a^\dagger a \right>",
		title = "Cavity occupation"
	)
	lines!(ax, Ndot_range, n_cav_ss_num; linestyle=:dash, alpha=0.7, label=L"DQD")
	lines!(ax, Ndot_range, n_cav_ss_ana; linestyle=:dash, alpha=0.7, label=L"No-DQD")
	axislegend(ax, position=:lt)
	if save_fig
	    save(joinpath(@__DIR__, "N_cav.png"), f)
	else
		display(f)
	end
end