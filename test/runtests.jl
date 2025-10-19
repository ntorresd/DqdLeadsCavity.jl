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
    dqd_leads = DqdLeads(Dqd(1. , 1., 1), Leads(2, 2, 2.))
    @test dqd_leads.dqd.blockade == true
    @test dqd_leads.dqd.U == Inf
    @test dqd_leads.leads.μ_avg == 0.0
end