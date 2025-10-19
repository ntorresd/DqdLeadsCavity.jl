mutable struct DqdLeads
    dqd::Dqd
    leads::Leads
end

function Ω(dqd_leads::DqdLeads)
    return Ω(dqd_leads.dqd)
end

function θ(dqd_leads::DqdLeads)
    return θ(dqd_leads.dqd)
end

function OnsiteEnergies(dqd_leads::DqdLeads)
    return OnsiteEnergies(dqd_leads.dqd)
end