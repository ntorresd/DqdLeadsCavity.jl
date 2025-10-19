mutable struct Cavity
    ωc::Real
end

function Base.show(io::IO, cavity::Cavity)
    print(io, "ωc=$(cavity.ωc)\n")
end