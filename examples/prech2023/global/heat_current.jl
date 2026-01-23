begin
	tc_list = logrange(1e-3, 1e1, 1000);
	dqdObj = deepcopy(dqd_leads)
	ρss_list = []
	for tc in tc_list
		dqdObj.dqd.tc = tc
		ρss = steadystate(
			build_H_dqd_ge(dqdObj.dqd),
			build_L_ops_dqd_gl(dqdObj)
		)
		push!(ρss_list, ρss)
	end
end

begin
    J_ss_ana_L = []
    J_ss_ana_R = []

    J_ss_num_L = []
    J_ss_num_R = []

    local i = 1
    for tc in tc_list
        dqdObj.dqd.tc = tc
        μL, μR = get_chemical_potentials(dqdObj.leads)
        H_dqd = build_H_dqd_LR(dqdObj.dqd)
        dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
        N_op = dL' * dL + dR' * dR
        L_ops = build_L_ops_dqd_gl(dqdObj)
        # analytical heat currents
        J_L, J_R = get_heat_current_gl(dqdObj)
        push!(J_ss_ana_L, J_L)
        push!(J_ss_ana_R, J_R)
        # numerical heat currents
        DLρ = D_sop(L_ops[1], ρss_list[i]) + D_sop(L_ops[2], ρss_list[i]) +
              D_sop(L_ops[3], ρss_list[i]) + D_sop(L_ops[4], ρss_list[i])
        DRρ = D_sop(L_ops[5], ρss_list[i]) + D_sop(L_ops[6], ρss_list[i]) +
              D_sop(L_ops[7], ρss_list[i]) + D_sop(L_ops[8], ρss_list[i])
        push!(J_ss_num_L, real(expect(H_dqd - μL * N_op, DLρ)))
        push!(J_ss_num_R, real(expect(H_dqd - μR * N_op, DRρ)))
        i += 1
    end
end

using Plots
using LaTeXStrings
Plots.pythonplot()
J_ss_plot_LR = Plots.plot(
    tc_list,
    [J_ss_ana_L, J_ss_ana_R],
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

J_ss_plot_L = Plots.plot(
    tc_list,
    [J_ss_ana_L, J_ss_num_L],
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

J_ss_plot_R = Plots.plot(
    tc_list,
    [J_ss_ana_R, J_ss_num_R],
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

plot(J_ss_plot_L, J_ss_plot_R, layout = (1, 2))