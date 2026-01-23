using Plots
using LaTeXStrings
Plots.pythonplot()

I_plot = Plots.plot(
    tc_list,
    [I_ana, I_num],
    xlabel = L"t_c",
    ylabel = L"I",
    label = [L"I^{ana}" L"I^{num}"],
    linestyle = [:solid :dash],
    linealpha = [0.3, 0.7],
    color = [:red :red],
    xaxis = :log,
    legend = :topright,
    dpi = 250
);