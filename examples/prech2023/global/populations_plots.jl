using Plots
# using PythonPlot
using LaTeXStrings
Plots.pythonplot()

n_dqd_ss_plot_ge = Plots.plot(
	tc_list,
	[n_dqd_ss_num_g, n_dqd_ss_ana_g, n_dqd_ss_num_e, n_dqd_ss_ana_e, α_ge_ss_num, α_ge_ss_ana],
	xlabel = L"t_c",
	ylabel = L"\left< d_\sigma^\dagger d_\sigma' \right>",
	label = [L"\bar{n}_g^{num}" L"\bar{n}_g^{ana}" L"\bar{n}_e^{num}" L"\bar{n}_e^{ana}" L"|α_{ge}^{num}|" L"|α_{ge}^{ana}|"],
	linestyle = [:dash :solid :dash :solid :dash :solid],
	linealpha = [0.7, 0.3, .7, 0.3, .7, 0.3],
	color = [:blue :blue :green :green :red :red],
	xaxis = :log,
	legend = :topright,
	dpi = 200
);

n_dqd_ss_plot_LR = Plots.plot(
	tc_list,
	[n_dqd_ss_num_L, n_dqd_ss_ana_L, n_dqd_ss_num_R, n_dqd_ss_ana_R, α_LR_ss_num, α_LR_ss_ana],
	xlabel = L"t_c",
	ylabel = L"\left< d_\alpha^\dagger d_\alpha \right>",
	label = [L"\bar{n}_L^{num}" L"\bar{n}_L^{ana}" L"\bar{n}_R^{num}" L"\bar{n}_R^{ana}" L"|α_{LR}^{num}|" L"|α_{LR}^{ana}|"],
	linestyle = [:dash :solid :dash :solid :dash :solid],
	linealpha = [0.7, 0.3, .7, 0.3, .7, 0.3],
	color = [:blue :blue :green :green :red :red],
	xaxis = :log,
	legend = :topright,
    dpi = 350
);

plot(n_dqd_ss_plot_ge, n_dqd_ss_plot_LR, layout = (1, 2))