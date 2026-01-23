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
	# ground state steady-state occupation
	n_dqd_ss_num_g = []
	n_dqd_ss_ana_g = []
	# excited state steady-state occupation
	n_dqd_ss_num_e = []
	n_dqd_ss_ana_e = []
	# ground-excited state coherence 
	α_ge_ss_num = []
	α_ge_ss_ana = []
	local i = 1
	for tc in tc_list
		dqdObj.dqd.tc = tc
		ρss = ρss_list[i]

		dg, de = build_dqd_fermi_ops_ge(dqdObj.dqd)
		ng, ne = dg' * dg, de' * de

		ng_ana, ne_ana = get_dqd_occupation_ge_gl(dqdObj)

		push!(n_dqd_ss_num_g, expect(ng, ρss))
		push!(n_dqd_ss_ana_g, ng_ana)
		push!(n_dqd_ss_num_e, expect(ne, ρss))
		push!(n_dqd_ss_ana_e, ne_ana)
		push!(α_ge_ss_num, abs(expect(de' * dg, ρss)))
		push!(α_ge_ss_ana, 0.0)
		i = i + 1
	end
end

begin
	# left dot steady-state occupation
	n_dqd_ss_num_L = []
	n_dqd_ss_ana_L = []
	# right dot steady-state occupation
	n_dqd_ss_num_R = []
	n_dqd_ss_ana_R = []
	# left-right state coherence 
	α_LR_ss_num = []
	α_LR_ss_ana = []
	local i = 1
	for tc in tc_list
		dqdObj.dqd.tc = tc
		ρss = ρss_list[i]

		dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
		nL, nR = dL' * dL, dR' * dR
		ket_0, ket_L, ket_R, ket_D = build_dqd_basis_LR(dqdObj.dqd)

		nL_ana, nR_ana = get_dqd_occupation_LR_gl(dqdObj)
		push!(n_dqd_ss_num_L, expect(nL, ρss))
		push!(n_dqd_ss_ana_L, nL_ana)
		push!(n_dqd_ss_num_R, expect(nR, ρss))
		push!(n_dqd_ss_ana_R, nR_ana)
		# push!(α_LR_ss_num, abs((expect(dR' * dL, ρss))))
		push!(α_LR_ss_num, abs(ket_L' * ρss * ket_R))
		push!(α_LR_ss_ana, get_dqd_coherence_LR_gl(dqdObj))
		i = i + 1
	end
end