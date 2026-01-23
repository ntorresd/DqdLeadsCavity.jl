using Plots
using LaTeXStrings
Plots.pythonplot()

J_plot_LR = Plots.plot(
    tc_list,
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
    tc_list,
    [J_ana_L_gl, J_num_L_gl],
    xlabel = L"t_c",
    ylabel = L"J_L",
    label = [L"J_L^{ana}" L"J_L^{num}"],
    linestyle = [:solid :dash],
    linealpha = [0.3, 0.7],
    color = [:blue :blue],
    xaxis = :log,
    legend = :topright,
    dpi = 250
);

J_plot_R = Plots.plot(
    tc_list,
    [J_ana_R_gl, J_num_R_gl],
    xlabel = L"t_c",
    ylabel = L"J_R",
    label = [L"J_R^{ana}" L"J_R^{num}"],
    linestyle = [:solid :dash],
    linealpha = [0.3, 0.7],
    color = [:red :red],
    xaxis = :log,
    legend = :topright,
    dpi = 250
);

plot(J_plot_L, J_plot_R, layout = (1, 2))