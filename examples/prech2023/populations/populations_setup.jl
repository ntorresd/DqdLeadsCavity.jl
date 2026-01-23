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

"""
Fermi factors for the global approach [Potts2021]
"""
function get_fermi_factors_gl(dqd_leads::DqdLeads)
	ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
	μL, μR = get_chemical_potentials(dqd_leads.leads)
	TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR

	fLg, fLe = fermi(ϵg, μL, TL), fermi(ϵe, μL, TL)
	fRg, fRe = fermi(ϵg, μR, TR), fermi(ϵe, μR, TR)

	return fLg, fLe, fRg, fRe
end

"""
Coupling strengths for the global approach according to [eq. (89) Potts2021]
"""
function get_coupling_strengths_gl(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	sθ2, cθ2 = sin(θ / 2.)^2, cos(θ / 2.)^2
	
	ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
	ΓLg, ΓLe = ΓL * sθ2, ΓL * cθ2
	ΓRg, ΓRe = ΓR * cθ2, ΓR * sθ2

	ΓLg, ΓLe, ΓRg, ΓRe
end

"""
Lindblad dissipator for the global approach according to [eq. (88) Potts2021]
"""
function build_L_ops_dqd_gl(dqd_leads::DqdLeads)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)

	dg, de = build_dqd_fermi_ops_ge(dqd_leads.dqd)
	
	L_ops = [
		sqrt(ΓLg * fLg) * dg', sqrt(ΓLg * (1. - fLg)) * dg,
		sqrt(ΓLe * fLe) * de', sqrt(ΓLe * (1. - fLe)) * de,
		sqrt(ΓRg * fRg) * dg', sqrt(ΓRg * (1. - fRg)) * dg,
		sqrt(ΓRe * fRe) * de', sqrt(ΓRe * (1. - fRe)) * de,
	]
	return L_ops
end
function get_L_ops_gl(dqd_leads_cavity::DqdLeadsCavityObj)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads_cavity.dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads_cavity.dqd_leads)

	dg, de = build_dqd_fermi_ops_ge(dqd_leads_cavity)
	
	L_ops = [
		sqrt(ΓLg * fLg) * dg', sqrt(ΓLg * (1. - fLg)) * dg,
		sqrt(ΓLe * fLe) * de', sqrt(ΓLe * (1. - fLe)) * de,
		sqrt(ΓRg * fRg) * dg', sqrt(ΓRg * (1. - fRg)) * dg,
		sqrt(ΓRe * fRe) * de', sqrt(ΓRe * (1. - fRe)) * de,
	]
	return L_ops
end

## Steady-state occupations (ground-excited)
"""
Analytical steady-state solution for the occupation of the DQD grond/excited state
according to the global approach[eq. (A27) Prech2023]
"""
function get_dqd_occupation_ge(dqd_leads::DqdLeads)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)

	ng = (ΓLg * fLg + ΓRg * fRg) / (ΓLg + ΓRg)
	ne = (ΓLe * fLe + ΓRe * fRe) / (ΓLe + ΓRe)

	return ng, ne
end


## Steady-state occupations (left-right)
"""
Analytical steady-state solution for the occupation of the DQD left and right states
according to global approach [eq. (A28) Prech2023]
"""
function get_dqd_occupation_LR(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	cθ2, sθ2 = cos(θ / 2.)^2, sin(θ / 2.)^2
	
	ng, ne = get_dqd_occupation_ge(dqd_leads)
	
	nL = cθ2 * ng + sθ2 * ne
	nR = sθ2 * ng + cθ2 * ne

	return nL, nR
end

"""
Analytical steady-state solution for the coherence between the ground and excited levels of the DQD according to the global approach [eq. (A29) Prech2023]
"""
function get_dqd_coherence_LR(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	ng, ne = get_dqd_occupation_ge(dqd_leads)
	α_abs = abs(sin(θ/2.) * cos(θ/2.) * (ng - ne))
	return(α_abs)
end


