export DqdLeadsCavityObj
export build_dqd_vladder_ops_LR, build_dqd_number_ops_LR
export build_dqd_σz_op
export build_cav_a_op

mutable struct DqdLeadsCavityObj
    dqd_leads::DqdLeads
    cavity::Cavity
end

function Base.show(io::IO, dqd_leads_cavity::DqdLeadsCavityObj)
    print(io,
        "Δϵ = $(dqd_leads_cavity.dqd_leads.dqd.Δϵ)\n",
        "ϵ_avg = $(dqd_leads_cavity.dqd_leads.dqd.ϵ_avg)\n",
        "tc = $(dqd_leads_cavity.dqd_leads.dqd.tc)\n",
        "γm = $(dqd_leads_cavity.dqd_leads.dqd.γm)\n",
        "γϕ = $(dqd_leads_cavity.dqd_leads.dqd.γϕ)\n",
        "U = $(dqd_leads_cavity.dqd_leads.dqd.U)\n",
        "Coulomb blockade: $(dqd_leads_cavity.dqd_leads.dqd.blockade)\n",
        "TL = $(dqd_leads_cavity.dqd_leads.leads.TL)\n",
        "TR = $(dqd_leads_cavity.dqd_leads.leads.TR)\n",
        "Δμ = $(dqd_leads_cavity.dqd_leads.leads.Δμ)\n",
        "μ_avg = $(dqd_leads_cavity.dqd_leads.leads.μ_avg)\n",
        "ΓL = $(dqd_leads_cavity.dqd_leads.ΓL)\n",
        "ΓR = $(dqd_leads_cavity.dqd_leads.ΓR)\n",
        "ωc = $(dqd_leads_cavity.cavity.ωc)\n",
        "κ = $(dqd_leads_cavity.cavity.κ)\n",
        "Tc = $(dqd_leads_cavity.cavity.Tc)\n",
        "ωd = $(dqd_leads_cavity.cavity.ωd)\n",
        "κin = $(dqd_leads_cavity.cavity.κin)\n",
        "Ndot = $(dqd_leads_cavity.cavity.Ndot)\n",
        "gJC = $(dqd_leads_cavity.cavity.gJC)\n",
        "Nc = $(dqd_leads_cavity.cavity.Nc)\n"
    )
end

# --- Derived parameters ---
@doc raw"""
DQD-driving detuning
"""
function get_Δd(dqd_leads_cavity::DqdLeadsCavityObj)
    Ω = get_Ω(dqd_leads_cavity.dqd_leads.dqd)
    return Ω - dqd_leads_cavity.cavity.ωd
end

# --- Operators ---
@doc raw"""
DQD vacuum ladder operators in the left-right basis
"""
function build_dqd_vladder_ops_LR(dqd::Dqd)
    if dqd.blockade
        ket_0, ket_L, ket_R = build_dqd_basis_LR(dqd)
        sL = ket_0 * ket_L'
        sR = ket_0 * ket_R'
    else
        sL = fdestroy(2, 2)
        sR = fdestroy(2, 1)
    end

    return sL, sR
end
function build_dqd_vladder_ops_LR(dqd_leads_cavity::DqdLeadsCavityObj)
    dqd = dqd_leads_cavity.dqd_leads.dqd
    sL_loc, sR_loc = build_dqd_vladder_ops_LR(dqd)
    id_cav = id_Cavity(dqd_leads_cavity.cavity)

    sL = tensor(sL_loc, id_cav)
    sR = tensor(sR_loc, id_cav)

    return sL, sR
end

@doc raw"""
Build DQD number operators in the LR basis
"""
function build_dqd_number_ops_LR(dqd::Dqd)
    if dqd.blockade
        ket_0, ket_L, ket_R = build_dqd_basis_LR(dqd)

        nL = ket_L * ket_L'
        nR = ket_R * ket_R'
        return nL, nR
    else
        nL = sL' * sL
        nR = sR' * sR
        nD = nL * nR
        return nL, nR, nD
    end
end
function build_dqd_number_ops_LR(dqd_leads_cavity::DqdLeadsCavityObj)
    dqd = dqd_leads_cavity.dqd_leads.dqd
    nL_loc, nR_loc = build_dqd_number_ops_LR(dqd)
    id_cav = id_Cavity(dqd_leads_cavity.cavity)

    nL = tensor(nL_loc, id_cav)
    nR = tensor(nR_loc, id_cav)

    return nL, nR
end

@doc raw"""
DQD vacuum ladder operators in the ground-excited basis
"""
function build_dqd_vladder_ops_ge(dqd::Dqd)
    sL, sR = build_dqd_vladder_ops_LR(dqd)
    θ = get_θ(dqd)

    sg = cos(θ/2.) * sL - sin(θ/2.) * sR # |0><g|
    se = sin(θ/2.) * sL + cos(θ/2.) * sR # |0><e|
    return sg, se
end
function build_dqd_vladder_ops_ge(dqd_leads_cavity::DqdLeadsCavityObj)
    sg_loc, se_loc = build_dqd_vladder_ops_ge(dqd_leads_cavity.dqd_leads.dqd)
    id_cav = id_Cavity(dqd_leads_cavity.cavity)

    sg = tensor(sg_loc, id_cav)
    se = tensor(se_loc, id_cav)

    return sg, se
end

@doc raw"""
DQD ground-excited ladder operators
"""
function build_dqd_ladder_ops_ge(dqd::Dqd)
	sg, se = build_dqd_vladder_ops_ge(dqd);

	σm = sg'*se;
	σp = se'*sg;
	return σm, σp
end
function build_dqd_ladder_ops_ge(dqd_leads_cavity::DqdLeadsCavityObj)
    sg, se = build_dqd_vladder_ops_ge(dqd_leads_cavity)

	σm = sg'*se;
	σp = se'*sg;
	return σm, σp   
end

@doc raw"""
DQD σz-operator in the ground-excited basis
"""
function build_dqd_σz_op(dqd::Dqd)
    sg, se = build_dqd_vladder_ops_ge(dqd);

    σz = se' * se - sg' * sg;
    return σz
end
function build_dqd_σz_op(dqd_leads_cavity::DqdLeadsCavityObj)
    sg, se = build_dqd_vladder_ops_ge(dqd_leads_cavity)

    σz = se' * se - sg' * sg;
    return σz
end

@doc raw"""
Cavity annihilation operator
"""
function build_cav_a_op(cavity::Cavity; disp::Bool = true)
    dim = cavity.Nc
    α = get_α(cavity)
    a = disp ? destroy(dim) + α : destroy(dim)
    return a
end
function build_cav_a_op(dqd_leads_cavity::DqdLeadsCavityObj; disp::Bool = true)
    a_loc = build_cav_a_op(dqd_leads_cavity.cavity; disp = disp)
    id_dqd = id_Dqd(dqd_leads_cavity.dqd_leads.dqd)

    a = tensor(id_dqd, a_loc)
    return a
end
