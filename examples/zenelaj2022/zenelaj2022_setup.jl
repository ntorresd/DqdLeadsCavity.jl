using Revise, DqdLeadsCavity

run_plots = true
save_fig = false

begin
	# DQD parameters
	local Δϵ = 17.57
	local ϵ_avg = 0.0
	local tc = 6.78
	local γm = 0.5 		# Relaxation rate
	local γphi = 3.92 	# Dephasing rates

	dqd = Dqd(Δϵ, ϵ_avg, tc, γm, γphi)

	# Leads
	local TL = 1e-7
	local TR = 1e-7
	local μL = 0.
	local μR = 0.

    # DQD-Leads
	local ΓL = 1.
	local ΓR = 1.

	leads = Leads(TL, TR, μL - μR, (μL + μR)/2.)
	dqd_leads = DqdLeads(dqd, leads, ΓL, ΓR)

	# Cavity parameters (from [zenelaj2022])
	local ωc = get_Ω(dqd)
	local κ = 0.337
	local Tc = 1e-7
	# Coherent driving
	local κin = 0.094
	local Ndot = 0.01
	local ωd = get_Ω(dqd)
	local gJC = 0.457

	cavity = Cavity(ωc, κ, Tc, ωd, κin, Ndot, gJC, 3)
	dqd_leads_cavity = DqdLeadsCavityObj(dqd_leads, cavity)
end

"""
Hamiltonian used in [zenelaj2022]
"""
function build_H_z2022(dqd_leads_cavity::DqdLeadsCavityObj)
	H_dqd_ge = build_H_dqd_ge(dqd_leads_cavity)
	H_drive_coh = build_H_cav_drive_coh(dqd_leads_cavity)
	H_JC = build_H_JC(dqd_leads_cavity)

	return H_dqd_ge + H_drive_coh + H_JC
end

"""
Lindblad operators for the DQD-cavity according to [zenelaj2022]
"""
function build_L_ops_z2022(dqd_leads_cavity::DqdLeadsCavityObj)
	# Derived parameters
	γm = dqd_leads_cavity.dqd_leads.dqd.γm
	γϕ = dqd_leads_cavity.dqd_leads.dqd.γϕ
	κ = dqd_leads_cavity.cavity.κ
	nth_cav = get_nth_cav(dqd_leads_cavity.cavity)

	# Derived operators
	σm, σp = build_dqd_ladder_ops_ge(dqd_leads_cavity)
	σz = build_dqd_σz_op(dqd_leads_cavity)
	a = build_cav_a_op(dqd_leads_cavity)

	# Dissipator
    L_ops_dqd = [
        build_L_ops_dqd_gl(dqd_leads_cavity);        
		[sqrt(γm) * σm, sqrt(γϕ)/2 * σz]        # relaxation and dephasing
    ]
    L_ops_cav = [		
        sqrt(κ * (nth_cav + 1)) * a, sqrt(κ * nth_cav) * a' # Cavity environment
    ]
    return([L_ops_dqd; L_ops_cav])
end

"""
Current in the large drive limit and the zero temperature limit
"""
function get_particle_current_ld(dqd_leads::DqdLeads)
    ΓLin, ΓRin, ΓLout, ΓRout = get_coupling_strengths_gl(dqd_leads);
    Γin = ΓLin + ΓRin;
    Γout = ΓLout + ΓRout;

    return (ΓLin * ΓRout - ΓRin * ΓLout)/(2 * Γin + Γout)
end

"""
Steady-state photton number in the cavity without DQD
"""
function get_n_cav_ss(cavity::Cavity)
	Δc = get_Δc(cavity)
	return abs(2. * sqrt(cavity.κin * cavity.Ndot)/(cavity.κ + 2im * Δc))^2
end

"""
Heat current in the large drive limit
"""
function get_heat_current_ld(dqd_leads_cavity::DqdLeadsCavityObj)
    ΓLin, ΓRin, ΓLout, ΓRout = get_coupling_strengths_gl(dqd_leads_cavity.dqd_leads);
    Γin = ΓLin + ΓRin;
    Γout = ΓLout + ΓRout;
	ϵg, ϵe = get_eigen_energies(dqd_leads_cavity.dqd_leads.dqd)

    return (ΓLin * Γout * ϵg - Γin * ΓLout * ϵe)/(2 * Γin + Γout)
end