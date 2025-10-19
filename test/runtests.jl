using DqdLeadsCavity
using Test

@testset "Test Dqd structure" begin
    @test Dqd(1., 1, 1., 1).blockade == false;
    @test Dqd(1, 1., 1).U == Inf;
    @test Dqd(1, 1, 1.).blockade == true;
end

@testset "Test Leads structure" begin
    @test Leads(1, 1., 1.).μ_avg == 0.0
end

@testset "Test DqdLeads structure" begin
    dqd = Dqd(1. , 1., 1)
    leads = Leads(2, 2, 2.)
    dqd_leads = DqdLeads(dqd, leads)

    @test dqd_leads.dqd.blockade == true
    @test dqd_leads.dqd.U == Inf
    @test dqd_leads.leads.μ_avg == 0.0

    @test Ω(dqd_leads.dqd) == Ω(dqd)
    @test θ(dqd_leads.dqd) == θ(dqd)
    @test OnsiteEnergies(dqd_leads.dqd) == OnsiteEnergies(dqd)
end