export Dqd
export get_Ω, get_θ, get_onsite_energies, get_eigen_energies
export build_dqd_basis_LR, build_dqd_ladder_ops_LR
export build_dqd_number_ops_LR, build_dqd_number_ops_ge, build_dqd_number_op
export build_dqd_basis_ge, build_dqd_ladder_ops_ge

# --- Dqd structure ---
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

# --- Derived parameters ---
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

# --- LR Basis ---

@doc raw"""
Build DQD L-R ket basis
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
Build DQD L-R ladder operators
"""
function build_dqd_ladder_ops_LR(dqd::Dqd)
    if dqd.blockade
        ket_0, ket_L, ket_R = build_dqd_basis_LR(dqd)
        cL = ket_0 * ket_L'
        cR = ket_0 * ket_R'
    else
        cL = fdestroy(2, 2)
        cR = fdestroy(2, 1)
    end

    return cL, cR
end

@doc raw"""
Build DQD number operators in the LR basis
"""
function build_dqd_number_ops_LR(dqd::Dqd)
    cL, cR = build_dqd_ladder_ops_LR(dqd)

    nL = cL' * cL
    nR = cR' * cR
    nD = nL * nR

    return nL, nR, nD
end

"""
Build DQD number operators in the g-e basis
"""
function build_dqd_number_ops_ge(dqd::Dqd)
    cg, ce = build_dqd_ladder_ops_ge(dqd)

    ng = cg' * cg
    ne = ce' * ce
    nD = ng * ne

    return ng, ne, nD
end

@doc raw"""
Build DQD total number operator
"""
function build_dqd_number_op(dqd::Dqd)
    nL, nR, nD = build_dqd_number_ops_LR(dqd)
    return nL + nR
end

# --- g-e basis ---

@doc raw"""
Build DQD g-e ket basis
"""
function build_dqd_basis_ge(dqd::Dqd)
    dqd.blockade || throw(
        ArgumentError(
            "DQD g-e basis representation is not supported for dqd.blockade=$(dqd.blockade)"
        ),
    )
    ket_0, ket_L, ket_R = build_dqd_basis_LR(dqd)
    θ = get_θ(dqd)

    ket_g = cos(θ/2) * ket_L - sin(θ/2) * ket_R
    ket_e = sin(θ/2) * ket_L + cos(θ/2) * ket_R

    return ket_0, ket_g, ket_e
end

@doc raw"""
Build DQD g-e ladder operators
"""
function build_dqd_ladder_ops_ge(dqd::Dqd)
    cL, cR = build_dqd_ladder_ops_LR(dqd)
    θ = get_θ(dqd)

    cg = cos(θ/2.) * cL - sin(θ/2.) * cR
    ce = sin(θ/2.) * cL + cos(θ/2.) * cR
    return cg, ce
end
