export DqdLeads
export get_transition_rates_ge, get_fermi_ge

mutable struct DqdLeads
    dqd::Dqd
    leads::Leads
    ΓL::Real
    ΓR::Real
end

function Base.show(io::IO, dqd_leads::DqdLeads)
    μL, μR = get_chemical_potentials(dqd_leads.leads)
    print(io,
        "Δϵ = $(dqd_leads.dqd.Δϵ)\n",
        "ϵ_avg = $(dqd_leads.dqd.ϵ_avg)\n",
        "tc = $(dqd_leads.dqd.tc)\n",
        "γm = $(dqd_leads.dqd.γm)\n",
        "γϕ = $(dqd_leads.dqd.γϕ)\n",
        "U = $(dqd_leads.dqd.U)\n",
        "Coulomb blockade: $(dqd_leads.dqd.blockade)\n",
        "TL = $(dqd_leads.leads.TL)\n",
        "TR = $(dqd_leads.leads.TR)\n",
        "Δμ = $(dqd_leads.leads.Δμ)\n",
        "μ_avg = $(dqd_leads.leads.μ_avg)\n",
        "μL = $μL\n",
        "μR = $μR\n",
        "ΓL = $(dqd_leads.ΓL)\n",
        "ΓR = $(dqd_leads.ΓR)\n"
    )
end

function get_onsite_energies(dqd_leads::DqdLeads)
    return get_onsite_energies(dqd_leads.dqd)
end

function get_chemical_potentials(dqd_leads::DqdLeads)
    return get_chemical_potentials(dqd_leads.leads)
end

function get_transition_rates_ge(dqd_leads::DqdLeads; side::Any = false)
    # Parameters
    ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
    μL, μR = get_chemical_potentials(dqd_leads.leads)
    TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR
    ΓL = dqd_leads.ΓL
    ΓR = dqd_leads.ΓR
    # Derived parameters
    θ = get_θ(dqd_leads.dqd)
    sθ2 = sin(θ/2)^2
    cθ2 = cos(θ/2)^2

	# Fermi functions (ground-excited)
	fLg = fermi(ϵg, μL, TL);
	fRg = fermi(ϵg, μR, TR);

	fLe = fermi(ϵe, μL, TL);
	fRe = fermi(ϵe, μR, TR);

	# Rates (ground-excited)
	Γg0L = ΓL * fLg * cθ2;
	Γg0R = ΓR * fRg * sθ2;

	Γe0L = ΓL * fLe * sθ2;
	Γe0R = ΓR * fRe * cθ2;

	Γ0gL = ΓL * (1. - fLg) * cθ2;
	Γ0gR = ΓR * (1. - fRg) * sθ2;

	Γ0eL = ΓL * (1. - fLe) * sθ2;
	Γ0eR = ΓR * (1. - fRe) * cθ2;

	if side == "left"
		return Γg0L, Γe0L, Γ0gL, Γ0eL
	elseif side == "right"
		return Γg0R, Γe0R, Γ0gR, Γ0eR
	else		
		Γg0 = Γg0L + Γg0R;
		Γe0 = Γe0L + Γe0R;
		Γ0g = Γ0gL + Γ0gR; 
		Γ0e = Γ0eL + Γ0eR;
		return Γg0, Γe0, Γ0g, Γ0e
	end
    # ΓLg, ΓLe, ΓRg, ΓRe (eq (89) [Potts et. al. 2021])
    # return ΓL * sθ2, ΓL * cθ2, ΓR * sθ2, ΓR * sθ2 
end

function get_fermi_ge(dqd_leads::DqdLeads)
    ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
    μL, μR = get_chemical_potentials(dqd_leads.leads)
    TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR

    # fLg, fLe, fRg, fRe
    return fermi(ϵg, μL, TL), fermi(ϵe, μL, TL), fermi(ϵg, μR, TR), fermi(ϵe, μR, TR)
end
