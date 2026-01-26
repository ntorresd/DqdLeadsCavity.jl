using Plots
using LaTeXStrings
Plots.pythonplot()

J_plot_LR = Plots.plot(
    tc_range,
    [J_ana_L_gl, J_ana_R_gl],
    xlabel = L"t_c",
    ylabel = L"J",
    label = [L"J_L^{ana}" L"J_R^{ana}"],
    linestyle = [:solid :solid],
    linealpha = [0.7, 0.7],
    color = [:blue :red],
    xaxis = :log,
    legend = :topright,
    dpi = 250
);

J_plot_L = Plots.plot(
    tc_range,
    [J_ana_L_gl, J_num_L_gl, J_L_neqgf],
    xlabel = L"t_c",
    ylabel = L"J_L",
    label = [L"J_L^{ana}" L"J_L^{num}" L"J_L^{neqgf}"],
    linestyle = [:solid :dash :solid],
    linealpha = [0.3, 0.7, 0.7],
    color = [:blue :blue :black],
    xaxis = :log,
    legend = :left,
    dpi = 250
);

J_plot_R = Plots.plot(
    tc_range,
    [J_ana_R_gl, J_num_R_gl, J_R_neqgf],
    xlabel = L"t_c",
    ylabel = L"J_R",
    label = [L"J_R^{ana}" L"J_R^{num}" L"J_R^{neqgf}"],
    linestyle = [:solid :dash :solid],
    linealpha = [0.3, 0.7, 0.7],
    color = [:red :red :black],
    xaxis = :log,
    legend = :left,
    dpi = 250
);

plot(J_plot_L, J_plot_R, layout = (1, 2))