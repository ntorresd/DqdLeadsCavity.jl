export Leads

mutable struct Leads
    TL::Real
    TR::Real
    Δμ::Real
    μ_avg::Real
end

function Leads(TL::Real, TR::Real, Δμ::Real)
    return Leads(TL, TR, Δμ, 0.0)
end

function Base.show(io::IO, leads::Leads)
    print(io,
        "TL = $(leads.TL)\n",
        "TR = $(leads.TR)\n",
        "Δμ = $(leads.Δμ)\n",
        "μ_avg = $(leads.μ_avg)\n"
    )
end
