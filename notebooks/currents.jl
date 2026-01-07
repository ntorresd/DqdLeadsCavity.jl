### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ ad0bcd24-af3b-11f0-8293-7f8392b1fa19
begin
    import Pkg
    # Activate the notebooks project (expects Project.toml in ./notebooks)
    Pkg.activate(@__DIR__)
    # Ensure the local package is available even if the Manifest is missing/outdated
    Pkg.develop(path=joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end

# ╔═╡ 8544f202-d781-4bea-9173-81e73ad0713a
# Load dependencies
begin
	using Revise
	using DqdLeadsCavity, QuantumToolbox
	using Plots, LaTeXStrings
end

# ╔═╡ 06f282ab-242c-44e1-baff-563d4a06093f
md"""
# Analytical currents
"""

# ╔═╡ bb44e975-905f-44d4-b34b-104619234ff9
md""" 
## Setup cells
"""

# ╔═╡ 7b8660d0-755d-4698-9bc4-704c3de96bc0
# Prech2023
begin
	# Leads
	Γ = 1.
	# T = 10. * Γ
	T = 10 * Γ
	TL = T
	TR = T
	# Δμ = 6.5 * T
	Δμ =  0.05
	
	# Dqd
	# Δϵ = 0.1 * Γ
	# ϵ_avg = 0.
	Δϵ = 0.
	ϵ_avg = 2. * Γ
	tc = 6.5 * Γ
end

# ╔═╡ 12d4a732-ca5a-4be7-b441-785a157183d3
begin
	dqd = Dqd(Δϵ, ϵ_avg, tc, false);
	leads = Leads(TL, TR, Δμ);
	dqd_leads = DqdLeads(dqd, leads, Γ, Γ)
	μL, μR = get_chemical_potentials(dqd_leads)
end

# ╔═╡ 7f601170-e008-4c27-887b-f7821238a16b
md"""
## Particle currents
"""

# ╔═╡ 64122fef-c1c8-4782-94ae-e64e2d6a19a6
begin
	N_range = 10000
	tc_range = logrange(1e-4, 25 * tc, N_range)
	# Compute NEQF current (non-interacting)
	I_neqfg = map(tc_range) do tc
	    dqd.tc = tc
	    current_particle_avg_neqgf(dqd_leads)
	end
	# Compute analytical current (local, non-interacting)
	I_thcl = map(tc_range) do tc
		dqd.tc = tc
		current_particle_avg_thcl(dqd_leads)
	end
end

# ╔═╡ d4da5b01-569f-4200-b6c6-71bf42c2a8f1
current_particle_plot = plot(
	tc_range / Γ, [I_neqfg, I_thcl],
    xlabel = L"t_c / \Gamma",
    ylabel = L"\langle I \rangle / \Gamma",
	label = [L"I_{neqgf}" L"I_{thcl}"],
	linestyle = :dash,
    legend = :topleft,
    xaxis = :log,
    dpi = 200
)

# ╔═╡ 220a4005-4585-4efa-8253-4b99ec5d7230
md"""
## Heat currents
"""

# ╔═╡ 82756943-c2be-4bfd-a123-c6d14b72d01d
begin
	# Compute NEQF current (non-interacting)
	J_neqfg = map(tc_range) do tc
	    dqd.tc = tc
	    current_heat_avg_neqgf(dqd_leads; left = false)
	end
	# Compute analytical heat current (local, non-interacting)
	J_thcl = map(tc_range) do tc
		dqd.tc = tc
		current_heat_avg_thcl(dqd_leads; left = false)
	end
	# Compute analytical heat current (global, non-interacting)
	J_thcg = map(tc_range) do tc
		dqd.tc = tc
		current_heat_avg_thcg_ana(dqd_leads; left = false)
	end
end

# ╔═╡ 46ac7fb7-dcca-40cf-a413-b8f85e33ce23
current_heat_plot = plot(
	tc_range / Γ, [J_neqfg, J_thcl, J_thcg],
    xlabel = L"t_c / \Gamma",
    ylabel = L"\langle J \rangle / \Gamma",
	label = [L"J_{neqgf}" L"J_{thcl}" L"J_{thcg}"],
	linestyle = :dash,
    legend = :left,
    xaxis = :log,
    dpi = 200
)

# ╔═╡ 11750161-0695-4db7-b99c-463ed663a1a5
savefig(current_heat_plot, "../../plots/currents/current_heat_ana_thcl_thcg_neqgf.png")

# ╔═╡ e81abe39-1a81-4b53-a713-580f5a784e4a
md"""
## Numerical currents
"""

# ╔═╡ 29df519b-7f48-4fdb-853d-2a46503596bd
function current_heat_tdbk(
    dqd_leads::DqdLeads,
    ρ::QuantumObject,
    H_tdbk::QuantumObject,
    L_op_α::QuantumObject,
    μ_α
)
    N_dqd = build_dqd_number_op(dqd_leads.dqd)
    return real(expect((H_tdbk - μ_α * N_dqd) * L_op_α, ρ))
end

# ╔═╡ 35b1d4cd-86b1-4fa9-b77e-70ca56e2afcc
begin
	J_thcg_num = map(tc_range) do tc
	    dqd_leads.dqd.tc = tc
	    H_dqd = build_H_dqd_ge(dqd_leads.dqd)
	    L_thcg = build_L_ops_thcg(dqd_leads)
	    ρss = steadystate(H_dqd, L_thcg)
	    current_heat_tdbk(dqd_leads, ρss, H_dqd, sum(L_thcg[1:2]), μL)
	end
end

# ╔═╡ 05ee0f43-f0f6-4243-8fb9-59dc2082383c
begin
	current_heat_plot_num = plot(
		tc_range / Γ, [J_neqfg, J_thcg, J_thcg_num],
	    xlabel = L"t_c / \Gamma",
	    ylabel = L"\langle I \rangle / \Gamma",
		label = [L"J_{neqgf}" L"J_{thcg}" L"J_{thcg}^{num}"],
		linestyle = [:dash :dash :solid],
		alpha = [1 1 0.7],
	    legend = :left,
	    xaxis = :log,
	    dpi = 200
	)
end

# ╔═╡ aee3f088-b708-4c5a-b8f6-aa4d386095d6
begin
	plot(
		tc_range / Γ, J_thcg_num,
	    xlabel = L"t_c / \Gamma",
	    ylabel = L"\langle I \rangle / \Gamma",
		label = L"J_{thcg}^{num}",
	    legend = :left,
	    xaxis = :log,
	    dpi = 200
	)
end

# ╔═╡ 15f07b68-4960-4309-ac82-153c944a8b0a
get_chemical_potentials(dqd_leads)

# ╔═╡ 4bd5c314-f02e-40d2-9442-77935240a0d3
context_LR(dqd_leads)

# ╔═╡ Cell order:
# ╟─06f282ab-242c-44e1-baff-563d4a06093f
# ╟─bb44e975-905f-44d4-b34b-104619234ff9
# ╠═ad0bcd24-af3b-11f0-8293-7f8392b1fa19
# ╠═8544f202-d781-4bea-9173-81e73ad0713a
# ╠═7b8660d0-755d-4698-9bc4-704c3de96bc0
# ╠═12d4a732-ca5a-4be7-b441-785a157183d3
# ╟─7f601170-e008-4c27-887b-f7821238a16b
# ╠═64122fef-c1c8-4782-94ae-e64e2d6a19a6
# ╠═d4da5b01-569f-4200-b6c6-71bf42c2a8f1
# ╟─220a4005-4585-4efa-8253-4b99ec5d7230
# ╠═82756943-c2be-4bfd-a123-c6d14b72d01d
# ╠═46ac7fb7-dcca-40cf-a413-b8f85e33ce23
# ╠═11750161-0695-4db7-b99c-463ed663a1a5
# ╟─e81abe39-1a81-4b53-a713-580f5a784e4a
# ╠═29df519b-7f48-4fdb-853d-2a46503596bd
# ╠═35b1d4cd-86b1-4fa9-b77e-70ca56e2afcc
# ╠═05ee0f43-f0f6-4243-8fb9-59dc2082383c
# ╠═aee3f088-b708-4c5a-b8f6-aa4d386095d6
# ╠═15f07b68-4960-4309-ac82-153c944a8b0a
# ╠═4bd5c314-f02e-40d2-9442-77935240a0d3
