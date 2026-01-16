export Cavity
export get_Δc, get_Δd, get_α
export get_nth_cav, get_n_cav_ss
export build_id_cavity

mutable struct Cavity
    ωc::Real    # Cavity resonance frequency
    κ::Real     # Cavity decay rate
    Tc::Real    # Cavity temperature
    ωd::Real    # Driving frequency
    κin::Real   # Driving coupling strength
    Ndot::Real  # Rate of impinging photons
    gJC::Real   # Jaynes-Cummings coupling strength
    Nc::Int     # Cavity Hilbert space dimension
end

function Cavity(ωc::Real, κ::Real, Tc::Real)
    return Cavity(ωc, κ, Tc, 0.0, 0.0, 0.0, 0.0, 10)
end

function Cavity(
    ωc::Real, κ::Real, Tc::Real,
    ωd::Real, κin::Real, Ndot::Real
)
    return Cavity(ωc, κ, Tc, ωd, κin, Ndot, 0.0, 10)
end

function Cavity(
    ωc::Real, κ::Real, Tc::Real,
    ωd::Real, κin::Real, Ndot::Real, gJC::Real
)
    return Cavity(ωc, κ, Tc, ωd, κin, Ndot, gJC, 10)
end

function Base.show(io::IO, cavity::Cavity)
    print(io, 
    "ωc = $(cavity.ωc)\n",
    "κ = $(cavity.κ)\n",
    "Tc = $(cavity.Tc)\n",
    "ωd = $(cavity.ωd)\n",
    "κin = $(cavity.κin)\n",
    "Ndot = $(cavity.Ndot)\n",
    "gJC = $(cavity.gJC)\n",
    "Nc = $(cavity.Nc)\n"
    )
end

# --- Derived parameters ---
@doc raw"""
DQD-cavity detuning
"""
function get_Δc(cavity::Cavity)
    return cavity.ωc - cavity.ωd
end

@doc raw"""
Cavity displacement field
"""
function get_α(cavity::Cavity)
	α = (-2im * sqrt(cavity.κin * cavity.Ndot)/(cavity.κ + 2im * get_Δc(cavity)));
	return α
end

@doc raw"""
Thermal occupation at the cavity frequency
"""
function get_nth_cav(cavity::Cavity)
    nth_cav = 1. / (exp(cavity.ωc/cavity.Tc) - 1.)
    return nth_cav
end

@doc raw"""
Steady-state photon number in the cavity without DQD
"""
function get_n_cav_ss(cavity::Cavity)
	n_cav_ss = abs(2. * sqrt(cavity.κin * cavity.Ndot)/(cavity.κ + 2im * get_Δc(cavity)))^2
	return n_cav_ss
end

# --- Cavity Operators ---
@doc raw"""
Cavity Identity operator
"""
function build_id_cavity(cavity::Cavity)
    dim = cavity.Nc
    return qeye(dim)
end
