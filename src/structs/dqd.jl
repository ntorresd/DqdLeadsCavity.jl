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