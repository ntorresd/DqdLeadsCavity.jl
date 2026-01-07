export DqdLeads
export get_transition_rates_ge, get_fermi_ge

mutable struct DqdLeads
    dqd::Dqd
    leads::Leads
    ΓL::Real
    ΓR::Real
end

function Base.show(io::IO, dqd_leads::DqdLeads)
    print(io,
        "Δϵ = $(dqd_leads.dqd.Δϵ)\n",
        "ϵ_avg = $(dqd_leads.dqd.ϵ_avg)\n",
        "tc = $(dqd_leads.dqd.tc)\n",
        "U = $(dqd_leads.dqd.U)\n",
        "Coulomb blockade: $(dqd_leads.dqd.blockade)\n",
        "TL = $(dqd_leads.leads.TL)\n",
        "TR = $(dqd_leads.leads.TR)\n",
        "Δμ = $(dqd_leads.leads.Δμ)\n",
        "μ_avg = $(dqd_leads.leads.μ_avg)\n",
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

function get_transition_rates_ge(dqd_leads::DqdLeads)
    θ = get_θ(dqd_leads.dqd)
    ΓL = dqd_leads.ΓL
    ΓR = dqd_leads.ΓR

    sθ2 = sin(θ/2)^2
    cθ2 = cos(θ/2)^2

    # ΓLg, ΓLe, ΓRg, ΓRe (eq (89) [Potts et. al. 2021])
    return ΓL * sθ2, ΓL * cθ2, ΓR * sθ2, ΓR * sθ2 
end

function get_fermi_ge(dqd_leads::DqdLeads)
    ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
    μL, μR = get_chemical_potentials(dqd_leads.leads)
    TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR

    # fLg, fLe, fRg, fRe
    return fermi(ϵg, μL, TL), fermi(ϵe, μL, TL), fermi(ϵg, μR, TR), fermi(ϵe, μR, TR)
end