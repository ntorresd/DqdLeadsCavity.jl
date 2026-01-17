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
Analytical steady-state solution for the occupation of the DQD grond/excited state according to global approach [eq. (A27) Prech2023]
"""
function get_dqd_occupation_ge(dqd_leads::DqdLeads)
	ΓLg, ΓLe, ΓRg, ΓRe = get_coupling_strengths_gl(dqd_leads)
	fLg, fLe, fRg, fRe = get_fermi_factors_gl(dqd_leads)

	ng = (ΓLg * fLg + ΓRg * fRg) / (ΓLg + ΓRg)
	ne = (ΓLe * fLe + ΓRe * fRe) / (ΓLe + ΓRe)

	return ng, ne
end

# ╔═╡ b546be59-30ed-4ce6-b217-0822ca85bafa
begin
	tc_list = logrange(1e-3, 1e1, 1000);
	dqdObj = deepcopy(dqd_leads)
	ρss_list = []
	for tc in tc_list
		dqdObj.dqd.tc = tc
		ρss = steadystate(
			build_H_dqd_ge(dqdObj.dqd),
			get_L_ops_dqd_gl(dqdObj)
		)
		push!(ρss_list, ρss)
	end
end

# ╔═╡ eb7a9560-735d-40b5-b387-eea392f3ca80
begin
	# ground state steady-state occupation
	n_dqd_ss_num_g = []
	n_dqd_ss_ana_g = []
	# excited state steady-state occupation
	n_dqd_ss_num_e = []
	n_dqd_ss_ana_e = []
	# ground-excited state coherence 
	α_ge_ss_num = []
	α_ge_ss_ana = []
	local i = 1
	for tc in tc_list
		dqdObj.dqd.tc = tc
		ρss = ρss_list[i]

		dg, de = build_dqd_fermi_ops_ge(dqdObj.dqd)
		ng, ne = dg' * dg, de' * de

		ng_ana, ne_ana = get_dqd_occupation_ge(dqdObj)

		push!(n_dqd_ss_num_g, expect(ng, ρss))
		push!(n_dqd_ss_ana_g, ng_ana)
		push!(n_dqd_ss_num_e, expect(ne, ρss))
		push!(n_dqd_ss_ana_e, ne_ana)
		push!(α_ge_ss_num, abs(expect(de' * dg, ρss)))
		push!(α_ge_ss_ana, 0.0)
		i = i + 1
	end
end

# ╔═╡ e01aedb5-f6f7-4a7e-af82-952b8a715725
n_dqd_ss_plot_ge = plot(
	tc_list,
	[n_dqd_ss_num_g, n_dqd_ss_ana_g, n_dqd_ss_num_e, n_dqd_ss_ana_e, α_ge_ss_num, α_ge_ss_ana],
	xlabel = L"t_c",
	ylabel = L"\left< d_\sigma^\dagger d_\sigma' \right>",
	label = [L"\bar{n}_g^{num}" L"\bar{n}_g^{ana}" L"\bar{n}_e^{num}" L"\bar{n}_e^{ana}" L"|α_{ge}^{num}|" L"|α_{ge}^{ana}|"],
	linestyle = [:dash :solid :dash :solid :dash :solid],
	linealpha = [.5, .1, .5, .1, .5, 0.1],
	color = [:blue :blue :green :green :red :red],
	xaxis = :log,
	legend = :right,
	dpi = 200
)

# ╔═╡ 9c37cdf9-f718-4e13-82d0-1c6906680e65
# savefig(n_dqd_ss_plot_ge, "../../plots/n_dqd_ss_plot_prech2023.png")

# ╔═╡ 4347ffe9-4619-43c2-a676-358ab5d6f17f
md"""
## Steady-state occupations (left-right)
"""

# ╔═╡ d50e14ac-ff1e-4309-a31b-b488d09dc359
"""
Analytical steady-state solution for the occupation of the DQD left and right states according to global approach [eq. (A28) Prech2023]
"""
function get_dqd_occupation_LR(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	cθ2, sθ2 = cos(θ)^2, sin(θ)^2
	
	ng, ne = get_dqd_occupation_ge(dqd_leads)
	
	nL = sθ2 * ng + cθ2 * ne
	nR = cθ2 * ng + sθ2 * ne

	return nL, nR
end

# ╔═╡ d6f22623-1b35-4a83-89eb-10686d758121
"""
Analytical steady-state solution for the coherence between the ground and excited levels of the DQD according to the global approach [eq. (A29) Prech2023]
"""
function get_dqd_coherence_LR(dqd_leads::DqdLeads)
	θ = get_θ(dqd_leads.dqd)
	ng, ne = get_dqd_occupation_ge(dqd_leads)
	α_abs = abs(sin(θ/2.) * cos(θ/2.) * (ne - ng))
	return(α_abs)
end

# ╔═╡ b0bc5fb9-1e57-498e-8d6e-32c2285cd70a
begin
	# ground state steady-state occupation
	n_dqd_ss_num_L = []
	n_dqd_ss_ana_L = []
	# excited state steady-state occupation
	n_dqd_ss_num_R = []
	n_dqd_ss_ana_R = []
	# ground-excited state coherence 
	α_LR_ss_num = []
	α_LR_ss_ana = []
	local i = 1
	for tc in tc_list
		dqdObj.dqd.tc = tc
		ρss = ρss_list[i]

		dL, dR = build_dqd_fermi_ops_LR(dqdObj.dqd)
		nL, nR = dL' * dL, dR' * dR
		ket_0, ket_L, ket_R, ket_D = build_dqd_basis_LR(dqdObj.dqd)

		nL_ana, nR_ana = get_dqd_occupation_LR(dqdObj)
		push!(n_dqd_ss_num_L, expect(nL, ρss))
		push!(n_dqd_ss_ana_L, nL_ana)
		push!(n_dqd_ss_num_R, expect(nR, ρss))
		push!(n_dqd_ss_ana_R, nR_ana)
		# push!(α_LR_ss_num, abs((expect(dR' * dL, ρss))))
		push!(α_LR_ss_num, abs(ket_L' * ρss * ket_R))
		push!(α_LR_ss_ana, get_dqd_coherence_LR(dqdObj))
		i = i + 1
	end
end

# ╔═╡ 00746111-b97f-4564-8767-6758dc1f3400
n_dqd_ss_plot_LR = plot(
	tc_list,
	[n_dqd_ss_num_L, n_dqd_ss_ana_L, n_dqd_ss_num_R, n_dqd_ss_ana_R, α_LR_ss_num, α_LR_ss_ana],
	xlabel = L"t_c",
	ylabel = L"\left< d_\sigma^\dagger d_\sigma' \right>",
	label = [L"\bar{n}_L^{num}" L"\bar{n}_L^{ana}" L"\bar{n}_R^{num}" L"\bar{n}_R^{ana}" L"|α_{LR}^{num}|" L"|α_{LR}^{ana}|"],
	linestyle = [:dash :solid :dash :solid :dash :solid],
	linealpha = [.5, .1, .5, .1, .5, 0.1],
	color = [:blue :blue :green :green :red :red],
	xaxis = :log,
	legend = :right,
	dpi = 200
)

# ╔═╡ 5264b7cb-0123-4d5d-b576-dc3a340c6573
plot(
	tc_list,
	[n_dqd_ss_num_L, n_dqd_ss_num_R, α_LR_ss_num],
	xlabel = L"t_c",
	ylabel = L"\left< d_\sigma^\dagger d_\sigma' \right>",
	label = [L"\bar{n}_L^{num}" L"\bar{n}_R^{num}" L"|α_{LR}^{num}|"],
	linestyle = [:dash :dash :dash],
	color = [:blue :green :red],
	xaxis = :log,
	legend = :right,
	dpi = 200
)

# ╔═╡ 982905b9-e7c2-40c5-9e8b-90800128bbd8
plot(
	tc_list,
	[n_dqd_ss_ana_L, n_dqd_ss_ana_R, α_LR_ss_ana],
	xlabel = L"t_c",
	ylabel = L"\left< d_\sigma^\dagger d_\sigma' \right>",
	label = [L"\bar{n}_L^{ana}" L"\bar{n}_R^{ana}" L"|α_{LR}^{ana}|"],
	linestyle = [:dash :dash :dash],
	color = [:blue :green :red],
	xaxis = :log,
	legend = :right,
	dpi = 200
)

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
# ╠═b546be59-30ed-4ce6-b217-0822ca85bafa
# ╠═eb7a9560-735d-40b5-b387-eea392f3ca80
# ╠═e01aedb5-f6f7-4a7e-af82-952b8a715725
# ╠═9c37cdf9-f718-4e13-82d0-1c6906680e65
# ╟─4347ffe9-4619-43c2-a676-358ab5d6f17f
# ╠═d50e14ac-ff1e-4309-a31b-b488d09dc359
# ╠═d6f22623-1b35-4a83-89eb-10686d758121
# ╠═b0bc5fb9-1e57-498e-8d6e-32c2285cd70a
# ╠═00746111-b97f-4564-8767-6758dc1f3400
# ╠═5264b7cb-0123-4d5d-b576-dc3a340c6573
# ╠═982905b9-e7c2-40c5-9e8b-90800128bbd8
