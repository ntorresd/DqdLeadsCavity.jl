export DqdLeads

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
        "ΓR = $(dqd_leads.ΓL)\n"
    )
end

function get_onsite_energies(dqd_leads::DqdLeads)
    return get_onsite_energies(dqd_leads.dqd)
end

function get_chemical_potentials(dqd_leads::DqdLeads)
    return get_chemical_potentials(dqd_leads.leads)
end