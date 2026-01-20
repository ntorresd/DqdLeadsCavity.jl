@testset "Test structures' defaults" begin
    # Dqd
    @test Dqd(1., 1, 1., 0., 0, 1.).blockade == false;
    @test Dqd(1, 1., 1, 0, 0.).U == Inf;
    @test Dqd(1, 1, 1., 0., 0.).blockade == true;

    # Leads
    @test Leads(1, 1., 1.).μ_avg == 0.0;

    # DqdLeads
    dqd = Dqd(1. , 1., 1, 0., 0.)
    leads = Leads(2, 2, 2.)
    dqd_leads = DqdLeads(dqd, leads, 1.0, 1.0)

    @test dqd_leads.dqd.blockade == true
    @test dqd_leads.dqd.U == Inf
    @test dqd_leads.leads.μ_avg == 0.0

end

@testset "Test parameters cross-structure consistency" begin
    dqd = Dqd(1. , 1., 1., 0., 0.)
    leads = Leads(2., 2., 2.)
    cavity = Cavity(0., 0., 0.)

    dqd_leads = DqdLeads(dqd, leads, 3.0, 3.0)
    dqd_leads_cavity = DqdLeadsCavityObj(dqd_leads, cavity)    

    @test get_Ω(dqd_leads.dqd) == get_Ω(dqd)
    @test get_θ(dqd_leads.dqd) == get_θ(dqd)
    @test get_Δd(dqd_leads_cavity) == get_Ω(dqd)
end