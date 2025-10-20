import Pkg
Pkg.activate(@__DIR__)
Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using DqdLeadsCavity

makedocs(
    sitename = "DqdLeadsCavity.jl",
    modules  = [DqdLeadsCavity],
    format   = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
    pages    = [
        "Home" => "index.md",
        "API" => "api.md"
    ],
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/ntorresd/DqdLeadsCavity.jl",
    devbranch = "main"
)