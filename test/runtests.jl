using QuantumToolbox
using DqdLeadsCavity
using Test

@testset "MyPackage.jl" begin
    include("test_structures.jl")
    include("test_sanity_checks.jl")
end
