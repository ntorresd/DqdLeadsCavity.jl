using QuantumToolbox

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
	n_dqd_num_g_gl = []
	n_dqd_ana_g_gl = []
	# excited state steady-state occupation
	n_dqd_num_e_gl = []
	n_dqd_ana_e_gl = []
	# ground-excited state coherence 
	α_ge_num_gl = []
	α_ge_ana_gl = []
	local i = 1
	for tc in tc_list
		dqdObj.dqd.tc = tc
		ρss = ρss_list[i]

		dg, de = build_dqd_fermi_ops_ge(dqdObj.dqd)
		ng, ne = dg' * dg, de' * de

		ng_ana, ne_ana = get_dqd_occupation_ge_gl(dqdObj)

		push!(n_dqd_num_g_gl, expect(ng, ρss))
		push!(n_dqd_ana_g_gl, ng_ana)
		push!(n_dqd_num_e_gl, expect(ne, ρss))
		push!(n_dqd_ana_e_gl, ne_ana)
		push!(α_ge_num_gl, abs(expect(de' * dg, ρss)))
		push!(α_ge_ana_gl, 0.0)
		i = i + 1
	end
end

begin
	# left dot steady-state occupation
	n_dqd_num_L_gl = []
	n_dqd_ana_L_gl = []
	# right dot steady-state occupation
	n_dqd_num_R_gl = []
	n_dqd_ana_R_gl = []
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
		push!(n_dqd_num_L_gl, expect(nL, ρss))
		push!(n_dqd_ana_L_gl, nL_ana)
		push!(n_dqd_num_R_gl, expect(nR, ρss))
		push!(n_dqd_ana_R_gl, nR_ana)
		# push!(α_LR_ss_num, abs((expect(dR' * dL, ρss))))
		push!(α_LR_ss_num, abs(ket_L' * ρss * ket_R))
		push!(α_LR_ss_ana, get_dqd_coherence_LR_gl(dqdObj))
		i = i + 1
	end
end