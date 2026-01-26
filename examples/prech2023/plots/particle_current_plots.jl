using Plots
using LaTeXStrings
Plots.pythonplot()

# I_plot = Plots.plot(
#     tc_range,
#     [I_L_ana_gl, I_L_num_gl],
#     xlabel = L"t_c",
#     ylabel = L"I",
#     label = [L"I^{gl}_{ana}" L"I^{gl}_{num}"],
#     linestyle = [:dashdot :dash :solid],
#     linealpha = [0.5, 0.7],
#     color = [:red :red],
#     xaxis = :log,
#     legend = :topright,
#     dpi = 250
# );

I_plot = Plots.plot(
    tc_range,
    [I_L_ana_gl, I_L_num_gl, I_L_neqgf],
    xlabel = L"t_c",
    ylabel = L"I",
    label = [L"I^{gl}_{ana}" L"I^{gl}_{num}" L"I^{NEGF}"],
    linestyle = [:dashdot :dash :solid],
    linealpha = [0.5, 0.7, 0.3],
    color = [:red :red :black],
    xaxis = :log,
    legend = :left,
    dpi = 250
);