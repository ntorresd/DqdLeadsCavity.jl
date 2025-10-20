export Cavity

mutable struct Cavity
    ωc::Real    # Cavity resonance frequency
    ωd::Real    # Driving frequency
end

function Cavity(ωc::Real)
    return Cavity(ωc, 0.0)
end

function Base.show(io::IO, cavity::Cavity)
    print(io, "ωc=$(cavity.ωc)\n")
end