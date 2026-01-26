using Plots
using LaTeXStrings
Plots.pythonplot()

I_plot_gl = Plots.plot(
    tc_range,
    [I_L_num_gl, I_L_ana_gl, I_L_neqgf],
    xlabel = L"t_c",
    ylabel = L"I_L",
    label = [L"I^{gl}_{num}" L"I^{gl}_{ana}" L"I^{NEGF}"],
    linestyle = [:dash :solid :solid],
    linealpha = [0.7, 0.5, 0.5],
    color = [:green :green :black],
    xaxis = :log,
    legend = :topright,
    dpi = 250
);

I_plot_loc = Plots.plot(
    tc_range,
    [I_L_ana_loc, I_L_num_loc, I_L_neqgf],
    xlabel = L"t_c",
    ylabel = L"I_L",
    label = [L"I^{loc}_{ana}" L"I^{loc}_{num}" L"I^{NEGF}"],
    linestyle = [:dash :solid :solid],
    linealpha = [0.7, 0.5, 0.5],
    color = [:green :green :black],
    xaxis = :log,
    legend = :topleft,
    dpi = 250
);

I_plot = Plots.plot(
    tc_range,
    [I_L_num_gl, I_L_num_loc, I_L_neqgf],
    xlabel = L"t_c",
    ylabel = L"I",
    label = [L"I^{gl}_{L}" L"I^{loc}_{L}" L"I^{NEGF}_{L}"],
    linestyle = [:dash :dash :solid],
    linealpha = [0.7, 0.7, 0.3],
    color = [:green :red :black],
    xaxis = :log,
    legend = :left,
    dpi = 250
);