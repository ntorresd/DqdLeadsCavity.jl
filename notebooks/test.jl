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
end

# ╔═╡ 7b8660d0-755d-4698-9bc4-704c3de96bc0
dqd = Dqd(1,1,1)

# ╔═╡ Cell order:
# ╟─ad0bcd24-af3b-11f0-8293-7f8392b1fa19
# ╠═8544f202-d781-4bea-9173-81e73ad0713a
# ╠═7b8660d0-755d-4698-9bc4-704c3de96bc0
