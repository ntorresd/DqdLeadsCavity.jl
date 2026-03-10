export get_fermi_factors_gl, get_coupling_strengths_gl
export build_L_ops_dqd_gl
export get_occupation_ge_gl_ana, get_occupation_LR_gl, get_coherence_LR_gl
export get_particle_current_gl
export get_heat_current_gl

@doc raw"""
Fermi factors for the global approach [Potts2021]
"""
function get_fermi_factors_gl(dqd_leads::DqdLeads)
	ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
	μL, μR = get_chemical_potentials(dqd_leads)
	TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR

	fLg, fLe = fermi(ϵg, μL, TL), fermi(ϵe, μL, TL)
	fRg, fRe = fermi(ϵg, μR, TR), fermi(ϵe, μR, TR)

	return fLg, fLe, fRg, fRe
end

@doc raw"""
Coupling strengths for the global approach according to
[eq. (A24) Prech2023] ≡ [eq. (89) Potts2021]
"""
function get_coupling_strengths_gl(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	cθ2, sθ2 = cos(θ / 2.)^2, sin(θ / 2.)^2
	
	ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
	# [eq. (A24) Prech2023] ≡ [eq. (89) Potts2021]
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
		sqrt(ΓLg * fLg) * dg', sqrt(ΓLg * (1. - fLg)) * dg, # L: |0> -> |g>, |g> -> |0>
		sqrt(ΓLe * fLe) * de', sqrt(ΓLe * (1. - fLe)) * de, # L: |0> -> |e>, |e> -> |0>
		sqrt(ΓRg * fRg) * dg', sqrt(ΓRg * (1. - fRg)) * dg, # R: |0> -> |g>, |g> -> |0>
		sqrt(ΓRe * fRe) * de', sqrt(ΓRe * (1. - fRe)) * de, # R: |0> -> |e>, |e> -> |0>
	]
	return L_ops
end

# Occupations and coherence
## Steady-state occupations (ground-excited)
@doc raw"""
Analytical steady-state solution for the occupation of the DQD grond/excited state
according to the global approach[eq. (A27) Prech2023]
"""
function get_occupation_ge_gl_ana(dqd_leads::DqdLeads)
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
function get_occupation_LR_gl(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	cθ2, sθ2 = cos(θ / 2.)^2, sin(θ / 2.)^2
	
	ng, ne = get_occupation_ge_gl_ana(dqd_leads)
	
	nL = cθ2 * ng + sθ2 * ne
	nR = sθ2 * ng + cθ2 * ne

	return nL, nR
end

@doc raw"""
Analytical steady-state solution for the coherence between the ground and excited levels of the DQD according to the global approach [eq. (A29) Prech2023]
"""
function get_coherence_LR_gl(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	ng, ne = get_occupation_ge_gl_ana(dqd_leads)
	α_abs = abs(sin(θ/2.) * cos(θ/2.) * (ng - ne))
	return(α_abs)
end

# Heat and particle currents
## Steady state particle currents
@doc raw"""
Analytical steady-state particle current for the non-interacting DQD
according to the global approach [eq. (A31) Prech2023].
Due to particle conservation I_R = -I_L.
"""
function get_particle_current_gl(dqd_leads::DqdLeads; left::Bool = true)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)
	_side = left ? 1. : -1.
	# particle currents
    I_L = (fLg - fRg) * ΓLg * ΓRg / (ΓLg + ΓRg) + (fLe - fRe) * ΓLe * ΓRe / (ΓLe + ΓRe)
	# I_R = (fRg - fLg) * ΓRg * ΓLg / (ΓRg + ΓLg) + (fRe - fLe) * ΓRe * ΓLe / (ΓRe + ΓLe)
	return _side * I_L
end

## Steady state heat currents
@doc raw"""
Analytical steady-state solution for the heat current
according to the global approach [eq. (A33) Prech2023] or
[eq. (B.8) Potts2021]
"""
function get_heat_current_gl(dqd_leads::DqdLeads; left::Bool = true)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)
	ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
	μL, μR = get_chemical_potentials(dqd_leads)

	Γg = ΓLg * ΓRg / (ΓLg + ΓRg)
	Γe = ΓLe * ΓRe / (ΓLe + ΓRe)

	Jg = left ? (ϵg - μL) * (fLg - fRg) * Γg : (ϵg - μR) * (fRg - fLg) * Γg
	Je = left ? (ϵe - μL) * (fLe - fRe) * Γe : (ϵe - μR) * (fRe - fLe) * Γe
	return Jg + Je
end
