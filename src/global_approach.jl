export get_fermi_factors_gl, get_coupling_strengths_gl
export build_L_ops_dqd_gl
export get_dqd_occupation_ge_gl, get_dqd_occupation_LR_gl, get_dqd_coherence_LR_gl
export get_particle_current_gl
export get_heat_current_gl

@doc raw"""
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

@doc raw"""
Coupling strengths for the global approach according to [eq. (89) Potts2021]
"""
function get_coupling_strengths_gl(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	cθ2, sθ2 = sin(θ / 2.)^2, cos(θ / 2.)^2
	
	ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
	ΓLg, ΓLe = ΓL * sθ2, ΓL * cθ2
	ΓRg, ΓRe = ΓR * cθ2, ΓR * sθ2

	ΓLg, ΓLe, ΓRg, ΓRe
end

@doc raw"""
Lindblad dissipator for the global approach according to [eq. (88) Potts2021].
This is valid for the non-interacting DQD case (U = 0).
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
function build_L_ops_dqd_gl(dqd_leads_cavity::DqdLeadsCavityObj)
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

# Occupations and coherence
## Steady-state occupations (ground-excited)
@doc raw"""
Analytical steady-state solution for the occupation of the DQD grond/excited state
according to the global approach[eq. (A27) Prech2023]
"""
function get_dqd_occupation_ge_gl(dqd_leads::DqdLeads)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)

	ng = (ΓLg * fLg + ΓRg * fRg) / (ΓLg + ΓRg)
	ne = (ΓLe * fLe + ΓRe * fRe) / (ΓLe + ΓRe)

	return ng, ne
end


## Steady-state occupations (left-right)
@doc raw"""
Analytical steady-state solution for the occupation of the DQD left and right states
according to global approach [eq. (A28) Prech2023]
"""
function get_dqd_occupation_LR_gl(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	cθ2, sθ2 = cos(θ / 2.)^2, sin(θ / 2.)^2
	
	ng, ne = get_dqd_occupation_ge_gl(dqd_leads)
	
	nL = cθ2 * ng + sθ2 * ne
	nR = sθ2 * ng + cθ2 * ne

	return nL, nR
end

@doc raw"""
Analytical steady-state solution for the coherence between the ground and excited levels of the DQD according to the global approach [eq. (A29) Prech2023]
"""
function get_dqd_coherence_LR_gl(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	ng, ne = get_dqd_occupation_ge_gl(dqd_leads)
	α_abs = abs(sin(θ/2.) * cos(θ/2.) * (ng - ne))
	return(α_abs)
end

# Heat and particle currents
## Steady state particle currents
@doc raw"""
Analytical steady-state solution for the particle current
according to the global approach [eq. (A31) Prech2023] 
"""
function get_particle_current_gl(dqd_leads::DqdLeads)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)

    I = (fLg - fRg) * ΓLg * ΓRg / (ΓLg + ΓRg) + (fLe - fRe) * ΓLe * ΓRe / (ΓLe + ΓRe)
end

## Steady state heat currents
@doc raw"""
Analytical steady-state solution for the heat current
according to the global approach [eq. (A33) Prech2023] or
[eq. (B.8) Potts2021]
"""
function get_heat_current_gl(dqd_leads::DqdLeads)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)
	ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
	μL, μR = get_chemical_potentials(dqd_leads.leads)

	Γg = ΓLg * ΓRg / (ΓLg + ΓRg)
	Γe = ΓLe * ΓRe / (ΓLe + ΓRe)

	JLg = (ϵg - μL) * (fLg - fRe) * Γg
	JLe = (ϵe - μL) * (fLe - fRe) * Γe

	JRg = (ϵg - μR) * (fRg - fLe) * Γg
	JRe = (ϵe - μR) * (fRe - fLe) * Γe

	return JLg + JLe, JRg + JRe
end
