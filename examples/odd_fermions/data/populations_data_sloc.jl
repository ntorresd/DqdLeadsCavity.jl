using QuantumToolbox

N_points = 10000
# Γ_range = logrange(1e-3, 50.0, N_points);
# Γ_range = logrange(1.0, 50.0, N_points);
# Γ_range = logrange(1e-3, 1e2, N_points);
Γ_range = logrange(0.5, 1e3, N_points)

begin
	local dqdObj = deepcopy(dqd_leads)
	# left dot steady-state occupation
	n0_num_loc = []
	n0_ana_loc = []
	n0_num_odd = []

	# left dot steady-state occupation
	nL_num_loc = []
	nL_ana_loc = []
	nL_num_odd = []
	# right dot steady-state occupation
	nR_num_loc = []
	nR_ana_loc = []
	nR_num_odd = []
	# left-right state coherence 
	αLR_num_loc = []
	αLR_ana_loc = []
	αLR_num_odd = []
	local i = 1
	for Γ in Γ_range
		# wrt left rate
		# dqdObj.ΓL = Γ
		# wrt rate proportion
		dqdObj.ΓL = dqdObj.ΓR * Γ
		
		local ket_0, ket_L, ket_R= build_dqd_basis_LR(dqdObj.dqd)
		local dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
		local n0 = ket_0 * ket_0'
		local nL, nR = dL' * dL, dR' * dR
		local H_dqd_LR = build_H_dqd_LR(dqdObj.dqd)
		local L_ops_loc = build_L_ops_dqd_loc(dqdObj)
		# Numeric populations local-blockade
		local ρss = steadystate(
			build_H_dqd_LR(dqdObj.dqd),
			build_L_ops_dqd_loc(dqdObj)
		)
		push!(n0_num_loc, real(expect(n0, ρss)))
		push!(nL_num_loc, real(expect(nL, ρss)))
		push!(nR_num_loc, real(expect(nR, ρss)))
		push!(αLR_num_loc, abs(ket_L' * ρss * ket_R))
		# Analytical populations semilocal-blockade
		n0_ana, nL_ana, nR_ana, nD_ana, αLR_ana = populations_ana_semilocal(dqdObj)
		push!(n0_ana_loc, n0_ana)
		push!(nL_ana_loc, nL_ana)
		push!(nR_ana_loc, nR_ana)
		push!(αLR_ana_loc, abs(αLR_ana))
		# Numerical populations odd-fermions
		local L_odd = liouvillian_odd(H_dqd_LR, L_ops_loc)
		local ρss_odd = steadystate(L_odd)
		push!(n0_num_odd, real(expect(n0, ρss_odd)))
		push!(nL_num_odd, real(expect(nL, ρss_odd)))
		push!(nR_num_odd, real(expect(nR, ρss_odd)))
		push!(αLR_num_odd, abs(ket_L' * ρss_odd * ket_R))
	end
end