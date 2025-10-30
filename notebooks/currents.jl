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
	using DqdLeadsCavity
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
	T = 10. * Γ
	Δμ = 6.5 * T
	
	# Dqd
	Δϵ = 0.1 * Γ
	ϵ_avg = 0.
	tc = 6.5
end

# ╔═╡ 12d4a732-ca5a-4be7-b441-785a157183d3
begin
	dqd = Dqd(Δϵ, ϵ_avg, tc);
	leads = Leads(T, T, Δμ);
	dqd_leads = DqdLeads(dqd, leads, Γ, Γ)
end

# ╔═╡ 7f601170-e008-4c27-887b-f7821238a16b
md"""
## Particle currents
"""

# ╔═╡ 64122fef-c1c8-4782-94ae-e64e2d6a19a6
begin
	N_range = 10000
	tc_range = logrange(0.01, 25 * tc, N_range)
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
	tc_range, [I_neqfg, I_thcl],
    xlabel = L"t_c / \Gamma",
    ylabel = L"\langle I \rangle / \Gamma",
	label = [L"I_{neqgf}" L"I_{thcl}"],
    legend = :topleft,
    xaxis = :log,
    dpi = 200
)

# ╔═╡ 964aa79a-1cdf-4bc6-9184-be657632f878
dqd.Δϵ

# ╔═╡ 9c43dfc4-30f4-4ed6-80a0-d7bda5251096
ϵL, ϵR = get_onsite_energies(dqd)

# ╔═╡ 220a4005-4585-4efa-8253-4b99ec5d7230
md"""
## Heat currents
"""

# ╔═╡ 82756943-c2be-4bfd-a123-c6d14b72d01d
begin
	# Compute NEQF current (non-interacting)
	J_neqfg = map(tc_range) do tc
	    dqd.tc = tc
	    current_heat_avg_neqgf(dqd_leads)
	end
	# Compute analytical current (local, non-interacting)
	J_thcl = map(tc_range) do tc
		dqd.tc = tc
		current_heat_avg_thcl(dqd_leads)
	end
end

# ╔═╡ 46ac7fb7-dcca-40cf-a413-b8f85e33ce23
current_heat_plot = plot(
	tc_range, [J_neqfg, J_thcl],
    xlabel = L"t_c / \Gamma",
    ylabel = L"\langle I \rangle / \Gamma",
	label = [L"J_{neqgf}" L"J_{thcl}"],
    legend = :bottomleft,
    xaxis = :log,
    dpi = 200
)

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
# ╠═964aa79a-1cdf-4bc6-9184-be657632f878
# ╠═9c43dfc4-30f4-4ed6-80a0-d7bda5251096
# ╟─220a4005-4585-4efa-8253-4b99ec5d7230
# ╠═82756943-c2be-4bfd-a123-c6d14b72d01d
# ╠═46ac7fb7-dcca-40cf-a413-b8f85e33ce23
