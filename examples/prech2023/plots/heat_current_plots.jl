using Plots
using LaTeXStrings
Plots.pythonplot()

J_plot_L_gl = Plots.plot(
    tc_range,
    [J_num_L_gl, J_ana_L_gl, J_L_neqgf],
    xlabel = L"t_c",
    ylabel = L"J_L",
    label = [L"J_L^{num}" L"J_L^{ana}" L"J_L^{neqgf}"],
    linestyle = [:dash :solid :solid],
    linealpha = [0.7, 0.3, 0.3],
    color = [:green :green :black],
    xaxis = :log,
    legend = :left,
    dpi = 250
);

J_plot_L_loc = Plots.plot(
    tc_range,
    [J_num_L_loc, J_ana_L_loc, J_L_neqgf],
    xlabel = L"t_c",
    ylabel = L"J_L",
    label = [L"J_L^{num}" L"J_L^{ana}" L"J_L^{neqgf}"],
    linestyle = [:dash :solid :solid],
    linealpha = [0.7, 0.3, 0.3],
    color = [:red :red :black],
    xaxis = :log,
    legend = :left,
    dpi = 250
);

J_plot_L = Plots.plot(
    tc_range,
    [J_num_L_gl, J_num_L_loc, J_L_neqgf],
    xlabel = L"t_c",
    ylabel = L"J",
    label = [L"J^{gl}_{L}" L"J^{loc}_{L}" L"J^{neqgf}_{L}"],
    linestyle = [:dash :dash :solid],
    linealpha = [0.7, 0.7, 0.3],
    color = [:green :red :black],
    xaxis = :log,
    legend = :left,
    dpi = 250
);

J_plot_R = Plots.plot(
    tc_range,
    [J_num_R_gl, J_num_R_loc, J_R_neqgf],
    xlabel = L"t_c",
    ylabel = L"J",
    label = [L"J^{gl}_{R}" L"J^{loc}_{R}" L"J^{neqgf}_{R}"],
    linestyle = [:dash :dash :solid],
    linealpha = [0.7, 0.7, 0.3],
    color = [:green :red :black],
    xaxis = :log,
    legend = :left,
    dpi = 250
);

plot(J_plot_L, J_plot_R, layout = (1, 2))