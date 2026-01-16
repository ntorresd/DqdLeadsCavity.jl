### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 8ccb7f3a-f245-11f0-a164-31a0d8f97e07
begin
    import Pkg
    # Activate the notebooks project (expects Project.toml in ./notebooks)
    Pkg.activate(@__DIR__)
    # Ensure the local package is available even if the Manifest is missing/outdated
    Pkg.develop(path=joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end

# ╔═╡ 0c8f46a8-ed2c-45b3-8999-b347f16158c5
begin
	using QuantumToolbox
	using Plots, LaTeXStrings
	using DqdLeadsCavity
end

# ╔═╡ 155e635f-9f3e-4018-b4c7-b7fe3bd15c2e
md"""
# Setup
"""

# ╔═╡ 6feccb92-7803-42a1-b4a7-06bbc9a5ea66
begin
	# Leads
	Γ = 1.
	T = 10. * Γ
	Δμ = 6.5 * T
	μ_avg = 0.0

	leads = Leads(T, T, Δμ, μ_avg)

	# DQD parameters
	Δϵ = 0.1 * Γ
	ϵ_avg = 0.0
	tc = 6.5
	U = 0.0
	
	γm = 0.			# Relaxation rate
	γphi = 3.92 	# Dephasing rate	
	
	dqd = Dqd(Δϵ, ϵ_avg, tc, γm, γphi, U)

	# Cavity parameters (from [zenelaj2022])
	ωc = get_Ω(dqd)
	κ = 0.337
	Tc = 1e-7
	# Coherent driving
	κin = 0.094
	Ndot = 0.01
	ωd = get_Ω(dqd)
	gJC = 0.457

	cavity = Cavity(ωc, κ, Tc, ωd, κin, Ndot, gJC, 3)
	
	dqd_leads = DqdLeads(dqd, leads, Γ, Γ)
	dqd_leads_cavity = DqdLeadsCavityObj(dqd_leads, cavity)
end

# ╔═╡ 0ab44a42-8094-42e8-a861-cb6a4eb1d927
"""
Fermi factors for the global approach [Potts2021]
"""
function get_fermi_factors_gl(dqd_leads::DqdLeads)
	ϵg, ϵe = get_eigen_energies(dqd_leads.dqd)
	μL, μR = get_chemical_potentials(dqd_leads.leads)
	TL, TR = dqd_leads.leads.TL, dqd_leads.leads.TR

	fLg, fLe = fermi(ϵg, μL, TL), fermi(ϵe, μL, TL)
	fRg, fRe = fermi(ϵg, μR, TR), fermi(ϵe, μR, TR)

	return fLg, fLe, fRg, fRe
end

# ╔═╡ ec3d5fe8-b23c-4401-a88f-0bdaa05c9ee4
"""
Coupling strengths for the global approach according to [eq. (89) Potts2021]
"""
function get_coupling_strengths_gl(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	sθ2, cθ2 = sin(θ)^2, cos(θ)^2
	
	ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
	ΓLg, ΓLe = ΓL * sθ2, ΓL * cθ2
	ΓRg, ΓRe = ΓR * cθ2, ΓR * sθ2

	ΓLg, ΓLe, ΓRg, ΓRe
end

# ╔═╡ bd44c07d-41bd-40a3-9e60-adfde3269600
begin
	"""
	Lindblad dissipator for the global approach according to [eq. (88) Potts2021]
	"""
	function get_L_ops_dqd_gl(dqd_leads::DqdLeads)
		ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
		fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)
	
		dg, de = build_dqd_fermi_ops_ge(dqd_leads.dqd)
		
		L_ops = [
			sqrt(ΓLg * fLg) * dg', sqrt(ΓLg * (1. - fLg)) * dg,
			sqrt(ΓLe * fLe) * de', sqrt(ΓLe * (1. - fLe)) * de,
			sqrt(ΓRg * fRg) * dg', sqrt(ΓRg * (1. - fRg)) * dg,
			sqrt(ΓRe * fRe) * de', sqrt(ΓRe * (1. - fRe)) * de,
		]
		return L_ops
	end
	function get_L_ops_gl(dqd_leads_cavity::DqdLeadsCavityObj)
		ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
		fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)
	
		dg, de = build_dqd_fermi_ops_ge(dqd_leads_cavity)
		
		L_ops = [
			sqrt(ΓLg * fLg) * dg', sqrt(ΓLg * (1. - fLg)) * dg,
			sqrt(ΓLe * fLe) * de', sqrt(ΓLe * (1. - fLe)) * de,
			sqrt(ΓRg * fRg) * dg', sqrt(ΓRg * (1. - fRg)) * dg,
			sqrt(ΓRe * fRe) * de', sqrt(ΓRe * (1. - fRe)) * de,
		]
		return L_ops
	end
end

# ╔═╡ 23587f58-575a-46fd-ae1f-b818fc208e7f
md"""
## Steady-state occupations (ground-excited)
"""

# ╔═╡ 78ae1607-0877-4a9b-a3a9-000cd8f8edd1
"""
Analytical steady-state solution for the occupation of the DQD grond/excited state in the global approach according to [eq. (A27) Prech2023]
"""
function get_dqd_occupation_ge(dqd_leads::DqdLeads)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)

	ng_ss = (ΓLg * fLg + ΓRg * fRg) / (ΓLg + ΓRg)
	ne_ss = (ΓLe * fLe + ΓRe * fRe) / (ΓLe + ΓRe)

	return ng_ss, ne_ss
end

# ╔═╡ eb7a9560-735d-40b5-b387-eea392f3ca80
let
	tc_list = logrange(1e-5, 100, 1000);
	dqdObj = deepcopy(dqd_leads)
	
	n_dqd_ss_num_g = []
	n_dqd_ss_num_e = []
	n_dqd_ss_ana_g = []
	n_dqd_ss_ana_e = []
	for tc in tc_list
		dqdObj.dqd.tc = tc
		ρss = steadystate(
			build_H_dqd_ge(dqdObj.dqd),
			get_L_ops_dqd_gl(dqdObj)
		)

		dg, de = build_dqd_fermi_ops_ge(dqdObj.dqd)
		ng, ne = dg' * dg, de' * de

		ng_ana, ne_ana = get_dqd_occupation_ge(dqdObj)

		push!(n_dqd_ss_num_g, expect(ng, ρss))
		push!(n_dqd_ss_num_e, expect(ne, ρss))
		push!(n_dqd_ss_ana_g, ng_ana)
		push!(n_dqd_ss_ana_e, ne_ana)
	end

	n_dqd_ss_plot = plot(
		tc_list,
		[n_dqd_ss_num_g, n_dqd_ss_num_e, n_dqd_ss_ana_g, n_dqd_ss_ana_e],
	    xlabel = L"t_c",
	    ylabel = L"\left< d_\sigma^\dagger d_\sigma \right>",
		label = [L"\bar{n}_g^{num}" L"\bar{n}_e^{num}" L"\bar{n}_g^{ana}" L"\bar{n}_e^{ana}"],
	    linestyle = [:solid :solid :dash :dash],
    	linealpha = [0.3, 0.3, 1.0, 1.0],
		xaxis = :log,
	    legend = :right,
	    dpi = 200
	)
end

# ╔═╡ 999cfdc7-f186-4ea4-8678-c7afdcd5c8b4
L_ops = get_L_ops_dqd_gl(dqd_leads)

# ╔═╡ c6035371-442d-47ea-9928-0b59f20eb7dc
H_dqd_ge = build_H_dqd_ge(dqd_leads.dqd)

# ╔═╡ 7dcc603f-7726-4884-bcae-ed860628a7a1
L_ops[1]

# ╔═╡ Cell order:
# ╠═8ccb7f3a-f245-11f0-a164-31a0d8f97e07
# ╠═0c8f46a8-ed2c-45b3-8999-b347f16158c5
# ╟─155e635f-9f3e-4018-b4c7-b7fe3bd15c2e
# ╠═6feccb92-7803-42a1-b4a7-06bbc9a5ea66
# ╠═0ab44a42-8094-42e8-a861-cb6a4eb1d927
# ╠═ec3d5fe8-b23c-4401-a88f-0bdaa05c9ee4
# ╠═bd44c07d-41bd-40a3-9e60-adfde3269600
# ╟─23587f58-575a-46fd-ae1f-b818fc208e7f
# ╠═78ae1607-0877-4a9b-a3a9-000cd8f8edd1
# ╠═eb7a9560-735d-40b5-b387-eea392f3ca80
# ╠═999cfdc7-f186-4ea4-8678-c7afdcd5c8b4
# ╠═c6035371-442d-47ea-9928-0b59f20eb7dc
# ╠═7dcc603f-7726-4884-bcae-ed860628a7a1
