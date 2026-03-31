using Revise, DqdLeadsCavity
using QuantumToolbox
using CairoMakie, LaTeXStrings
using ColorSchemes

run_plots = true
save_fig = false

# Setup
begin
	# DQD parameters
	local Δϵ = 0.0
	local ϵ_avg = 0.0
	local tc = 6.78
	# local γm = 0.5 		# Relaxation rate
	# local γϕ = 3.92 	# Dephasing rates
	local γm = 0.0 		# Relaxation rate
	local γϕ = 0.0 	# Dephasing rates

	# dqd = Dqd(Δϵ, ϵ_avg, tc, γm, γϕ) 		# Coulomb-blockade
	dqd = Dqd(0.0, ϵ_avg, tc, γm, γϕ, 0.0)	# Non-interacting 

	# Leads
	# local T = 1e-2 # don't erase
	local T = 10.
	local TL = T
	local TR = T
	local eV = 100.
	local μ_avg = 0.0

    # DQD-Leads
	local Γ = 1.
	local ΓL = Γ
	local ΓR = Γ
	# local ΓL = 2.
	# local ΓR = 1.

	leads = Leads(TL, TR, eV, μ_avg)
	dqd_leads = DqdLeads(dqd, leads, ΓL, ΓR)
end
