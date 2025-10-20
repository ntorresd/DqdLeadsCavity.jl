using DqdLeadsCavity
using Test

@testset "Test structures' defaults" begin
    # Dqd
    @test Dqd(1., 1, 1., 1).blockade == false;
    @test Dqd(1, 1., 1).U == Inf;
    @test Dqd(1, 1, 1.).blockade == true;

    # Leads
    @test Leads(1, 1., 1.).μ_avg == 0.0;

    # DqdLeads
    dqd = Dqd(1. , 1., 1)
    leads = Leads(2, 2, 2.)
    dqd_leads = DqdLeads(dqd, leads, 1.0, 1.0)

    @test dqd_leads.dqd.blockade == true
    @test dqd_leads.dqd.U == Inf
    @test dqd_leads.leads.μ_avg == 0.0

    @test get_Ω(dqd_leads.dqd) == get_Ω(dqd)
    @test get_θ(dqd_leads.dqd) == get_θ(dqd)
end

@testset "Test parameters cross-structure consistency" begin
    dqd = Dqd(1. , 1., 1.)
    leads = Leads(2., 2., 2.)
    dqd_leads = DqdLeads(dqd, leads, 3.0, 3.0)
    cavity = Cavity(0., 0.)

    @test get_Ω(dqd_leads.dqd) == get_Ω(dqd)
    @test get_θ(dqd_leads.dqd) == get_θ(dqd)
    @test get_Δd(dqd, cavity) == get_Ω(dqd)
end
