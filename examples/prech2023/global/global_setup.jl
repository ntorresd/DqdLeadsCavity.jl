using QuantumToolbox
using Revise, DqdLeadsCavity

# Setup
begin
	# Leads
	Γ = 1.
	T = 10. * Γ
	Δμ = 6.5 * T
	μ_avg = 0.0

	leads = Leads(T, T, Δμ, μ_avg)

	# DQD parameters
	Δϵ = 0.1 * Γ
	ϵ_avg = 0.0
	tc = 6.5
	U = 0.0
	
	γm = 0.			# Relaxation rate
	γphi = 3.92 	# Dephasing rate	
	
	dqd = Dqd(Δϵ, ϵ_avg, tc, γm, γphi, U)

	# Cavity parameters (from [zenelaj2022])
	ωc = get_Ω(dqd)
	κ = 0.337
	Tc = 1e-7
	# Coherent driving
	κin = 0.094
	Ndot = 0.01
	ωd = get_Ω(dqd)
	gJC = 0.457

	cavity = Cavity(ωc, κ, Tc, ωd, κin, Ndot, gJC, 3)
	
	dqd_leads = DqdLeads(dqd, leads, Γ, Γ)
	dqd_leads_cavity = DqdLeadsCavityObj(dqd_leads, cavity)
end
