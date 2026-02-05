using QuantumToolbox

N_points = 10000
tc_range = logrange(0.01, 1e3, N_points);

begin
	local dqdObj = deepcopy(dqd_leads)
	# ground state steady-state occupation
	n_dqd_num_g_gl = []
	n_dqd_ana_g_gl = []
	# excited state steady-state occupation
	n_dqd_num_e_gl = []
	n_dqd_ana_e_gl = []
	# ground-excited state coherence 
	α_ge_num_gl = []
	α_ge_ana_gl = []
	for tc in tc_range
		dqdObj.dqd.tc = tc
		local ρss = steadystate(
			build_H_dqd_ge(dqdObj.dqd),
			build_L_ops_dqd_gl(dqdObj)
		)
		local dg, de = build_dqd_fermi_ops_ge(dqdObj.dqd)
		local ng, ne = dg' * dg, de' * de

		local ng_ana, ne_ana = get_occupation_ge_gl(dqdObj)

		push!(n_dqd_num_g_gl, real(expect(ng, ρss)))
		push!(n_dqd_ana_g_gl, ng_ana)
		push!(n_dqd_num_e_gl, real(expect(ne, ρss)))
		push!(n_dqd_ana_e_gl, ne_ana)
		push!(α_ge_num_gl, abs(expect(de' * dg, ρss)))
		push!(α_ge_ana_gl, 0.0)
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
	α_LR_num_gl = []
	α_LR_ana_gl = []
	local i = 1
	for tc in tc_range
		dqdObj.dqd.tc = tc
		local ρss = steadystate(
			build_H_dqd_ge(dqdObj.dqd),
			build_L_ops_dqd_gl(dqdObj)
		)

		local dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
		local nL, nR = dL' * dL, dR' * dR
		local ket_0, ket_L, ket_R, ket_D = build_dqd_basis_LR(dqdObj.dqd)

		nL_ana, nR_ana = get_occupation_LR_gl(dqdObj)
		push!(n_dqd_num_L_gl, real(expect(nL, ρss)))
		push!(n_dqd_ana_L_gl, nL_ana)
		push!(n_dqd_num_R_gl, real(expect(nR, ρss)))
		push!(n_dqd_ana_R_gl, nR_ana)
		# push!(α_LR_ss_num, abs((expect(dR' * dL, ρss))))
		push!(α_LR_num_gl, abs(ket_L' * ρss * ket_R))
		push!(α_LR_ana_gl, get_coherence_LR_gl(dqdObj))
	end
end