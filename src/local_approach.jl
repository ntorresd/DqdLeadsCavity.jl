export get_fermi_factors_loc
export build_L_ops_dqd_loc
export get_particle_current_loc
export get_heat_current_loc

@doc raw"""
Fermi factors for the local approach [Potts2021]
"""
function get_fermi_factors_loc(dqd_leads::DqdLeads)
    ϵ_avg = dqd_leads.dqd.ϵ_avg
    μL, μR = get_chemical_potentials(dqd_leads)
	TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR

    fL, fR = fermi(ϵ_avg, μL, TL), fermi(ϵ_avg, μR, TR)
	return fL, fR
end

@doc raw"""
Lindblad dissipator for the non-interacting (U = 0) DQD
according to the local approach [eq. (92) Potts2021].
"""
function build_L_ops_dqd_loc(dqd_leads::DqdLeads)
    ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
    fL, fR = get_fermi_factors_loc(dqd_leads)
	dL, dR = build_dqd_fermi_ops_LR(dqd_leads.dqd)
	
	L_ops = [
		sqrt(ΓL * fL) * dL', sqrt(ΓL * (1. - fL)) * dL,
		sqrt(ΓR * fR) * dR', sqrt(ΓR * (1. - fR)) * dR
	]
	return L_ops
end
function build_L_ops_dqd_loc(dqd_leads_cavity::DqdLeadsCavityObj)
    ΓL, ΓR = dqd_leads_cavity.dqd_leads.ΓL, dqd_leads_cavity.dqd_leads.ΓR
    fL, fR = get_fermi_factors_loc(dqd_leads_cavity.dqd_leads)
	dL, dR = build_dqd_fermi_ops_LR(dqd_leads_cavity)
	
	L_ops = [
		sqrt(ΓL * fL) * dL', sqrt(ΓL * (1. - fL)) * dL,
		sqrt(ΓR * fR) * dR', sqrt(ΓR * (1. - fR)) * dR
	]
	return L_ops
end

@doc raw"""
Analytical steady-state particle current for the non-interacting DQD
according to the local approach[eq. (A12) Prech2023].
Due to conservation of particles I_R = -I_L
"""
function get_particle_current_loc(dqd_leads::DqdLeads; left::Bool = true)
    # Parameters
    tc = dqd_leads.dqd.tc
    ϵ_avg = dqd_leads.dqd.ϵ_avg
    Δϵ = dqd_leads.dqd.Δϵ
    μL, μR = get_chemical_potentials(dqd_leads)
    TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR
    ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
    # Fermi factors
    fL = fermi(ϵ_avg, μL, TL)
    fR = fermi(ϵ_avg, μR, TR)
    # current
    _side = left ? 1. : -1.
    numerator = 4 * tc^2 * (ΓL * ΓR) * (ΓL + ΓR)
    denominator = (ΓL + ΓR)^2 * (4. * tc^2 + (ΓL * ΓR) + 4. * (ΓL * ΓR) * Δϵ^2)
    return _side * (fL - fR) * numerator / denominator
end

@doc raw"""
Aalytical steady-state heat current through the non-interacting DQD
according to the local approach [eq. (B.6) Potts et. al. 2021]
"""
function get_heat_current_loc(dqd_leads::DqdLeads; left::Bool = true)
    # Parameters
    ϵ_avg = dqd_leads.dqd.ϵ_avg
    Δϵ = dqd_leads.dqd.Δϵ
    tc = dqd_leads.dqd.tc
    ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
    μL, μR = get_chemical_potentials(dqd_leads)
    TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR
    # Fermi factors
    fL = fermi(ϵ_avg, μL, TL)
    fR = fermi(ϵ_avg, μR, TR)
    factor = left ? (ϵ_avg - μL) * (fL - fR) : (ϵ_avg - μR) * (fR - fL)
    # effective rates
    Γ_avg = (ΓL + ΓR) / 2.
    γ = sqrt(ΓL * ΓR / Γ_avg^2) # Adimensional constant

    numerator = 2 * Γ_avg * (γ * tc)^2
    denominator = 4 * tc^2 + (γ * Γ_avg)^2 + (γ * Δϵ)^2
    return factor * numerator / denominator
end
