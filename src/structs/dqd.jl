mutable struct Dqd
    Δϵ::Real
    ϵ_avg::Real
    tc::Real
    U::Real
    blockade::Bool
end

function Dqd(Δϵ::Real, ϵ_avg::Real, tc::Real, U::Real)
    return Dqd(Δϵ, ϵ_avg, tc, U, false)
end

function Dqd(Δϵ::Real, ϵ_avg::Real, tc::Real)
    return Dqd(Δϵ, ϵ_avg, tc, Inf, true)
end

function Base.show(io::IO, dqd::Dqd)
    print(io,
        "Δϵ = $(dqd.Δϵ)\n",
        "ϵ_avg = $(dqd.ϵ_avg)\n",
        "tc = $(dqd.tc)\n",
        "U = $(dqd.U)\n",
        "Coulomb blockade: $(dqd.blockade)\n"
    )
end

function Ω(dqd::Dqd)
    return sqrt(4 * dqd.tc^2 + dqd.Δϵ^2)
end

function θ(dqd::Dqd)
    return atan(2. * dqd.tc / dqd.Δϵ) / 2.
end

function OnsiteEnergies(dqd::Dqd)
    ϵL = dqd.ϵ_avg + dqd.Δϵ / 2.
    ϵR = dqd.ϵ_avg - dqd.Δϵ / 2.
    return ϵL, ϵR
end