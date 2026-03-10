using Revise, DqdLeadsCavity
using QuantumToolbox

run_plots = true
save_fig = false

# Setup
begin
	# Leads
	local Γ = 1.
	local T = 100. * Γ
	# local Δμ = 6.5 * T
	local Δμ = 1e3
	local μ_avg = 0.0

	leads = Leads(T, T, Δμ, μ_avg)

	# DQD parameters
	# local Δϵ = 0. * Γ
	local Δϵ = 131.29 # Max concurrence for eV=650
	local ϵ_avg = 0.0
	# local tc = 0.0
	local tc = T / 2.
	local U = 0.0
	
	local γm = 0.			# Relaxation rate
	local γphi = 3.92 	# Dephasing rate

	dqd = Dqd(Δϵ, ϵ_avg, tc, γm, γphi, U)

	# Cavity parameters (from [zenelaj2022])
	local ωc = get_Ω(dqd)
	local κ = 0.337
	local Tc = 1e-7
	# Coherent driving
	local κin = 0.094
	local Ndot = 0.01
	local ωd = get_Ω(dqd)
	local gJC = 0.457

	cavity = Cavity(ωc, κ, Tc, ωd, κin, Ndot, gJC, 3)
	
	dqd_leads = DqdLeads(dqd, leads, Γ, Γ)
	dqd_leads_cavity = DqdLeadsCavityObj(dqd_leads, cavity)
end
