export Dqd
export get_dim, get_Ω, get_θ, get_onsite_energies, get_eigen_energies
export id_Dqd
export build_dqd_basis_LR
export build_dqd_vladder_ops_ge, build_dqd_ladder_ops_ge, build_dqd_σz_ge
export build_dqd_number_ops_ge, build_dqd_number_op

# --- Dqd structure ---
mutable struct Dqd
    Δϵ::Real
    ϵ_avg::Real
    tc::Real
    γm::Real
    γϕ::Real
    U::Real
    blockade::Bool
end

function Dqd(Δϵ::Real, ϵ_avg::Real, tc::Real, γm::Real, γϕ::Real, U::Real)
    return Dqd(Δϵ, ϵ_avg, tc, γm, γϕ, U, false)
end

function Dqd(Δϵ::Real, ϵ_avg::Real, tc::Real, γm::Real, γϕ::Real)
    return Dqd(Δϵ, ϵ_avg, tc, γm::Real, γϕ::Real, Inf, true)
end

function Base.show(io::IO, dqd::Dqd)
    print(io,
        "Δϵ = $(dqd.Δϵ)\n",
        "ϵ_avg = $(dqd.ϵ_avg)\n",
        "tc = $(dqd.tc)\n",
        "γm = $(dqd.γm)\n",
        "γϕ = $(dqd.γϕ)\n",
        "U = $(dqd.U)\n",
        "Coulomb blockade: $(dqd.blockade)\n"
    )
end

# --- Derived parameters ---
@doc raw"""
Get DQD Hilbert space dimension
"""
function get_dim(dqd::Dqd)
    dim = dqd.blockade ? 3 : 4
    return dim
end

@doc raw"""
DQD energy split in the g-e basis
"""
function get_Ω(dqd::Dqd)
    return sqrt(4 * dqd.tc^2 + dqd.Δϵ^2)
end

@doc raw"""
DQD mixing angle in the g-e basis
"""
function get_θ(dqd::Dqd)
    θ = acos(- dqd.Δϵ / get_Ω(dqd))
    return θ
end

@doc raw"""
DQD on-site energies
"""
function get_onsite_energies(dqd::Dqd)
    ϵL = dqd.ϵ_avg + dqd.Δϵ / 2.
    ϵR = dqd.ϵ_avg - dqd.Δϵ / 2.
    return ϵL, ϵR
end

@doc raw"""
DQD eigen-energies
"""
function get_eigen_energies(dqd::Dqd)
    Ω = get_Ω(dqd)
    ϵg = dqd.ϵ_avg - Ω / 2.
    ϵe = dqd.ϵ_avg + Ω / 2.
    return ϵg, ϵe
end

# --- DQD Operators ---
@doc raw"""
DQD Identity operator
"""
function id_Dqd(dqd::Dqd)
    dim = get_dim(dqd)
    return qeye(dim)
end

@doc raw"""
DQD L-R ket basis
"""
function build_dqd_basis_LR(dqd::Dqd)
    dim = get_dim(dqd)
    if dqd.blockade
        # DQD basis
        ket_0 = basis(dim, 0)
        ket_L = basis(dim, 1)
        ket_R = basis(dim, 2)

        return ket_0, ket_L, ket_R
    else
        # Dimensions
        dims = (2, 2)

        # DQD basis
        ket_0 = fock(dim, 0, dims = dims)  # |0,0⟩
        ket_L = fock(dim, 1, dims = dims)   # |0,1⟩
        ket_R = fock(dim, 2, dims = dims)   # |1,0⟩
        ket_D = fock(dim, 3, dims = dims)     # |1,1⟩

        return ket_0, ket_L, ket_R, ket_D
    end
end

@doc raw"""
DQD-σz operator in the ground-excited basis
"""
function build_dqd_σz_ge(dqd::Dqd)
	sg, se = build_dqd_vladder_ops_ge(dqd);
	
	σz = se'*se - sg'*sg;
	return σz
end

"""
DQD number operators in the g-e basis
"""
function build_dqd_number_ops_ge(dqd::Dqd)
    cg, ce = build_dqd_ladder_ops_ge(dqd)

    ng = cg' * cg
    ne = ce' * ce
    nD = ng * ne

    return ng, ne, nD
end

@doc raw"""
DQD total number operator
"""
function build_dqd_number_op(dqd::Dqd)
    nL, nR, nD = build_dqd_number_ops_LR(dqd)
    return nL + nR
end
