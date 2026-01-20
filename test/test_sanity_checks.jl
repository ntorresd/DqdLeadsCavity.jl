@testset "Test DQD basis and operators in the left-right basis (no-blockade)" begin
    dqd = Dqd(17.57, 0.0, 6.78, 0., 0., 0.)
    ket_0, ket_L, ket_R, ket_LR = build_dqd_basis_LR(dqd)
    dL, dR = build_dqd_fermi_ops_LR(dqd)

    @test dL * ket_L ≈ ket_0
    @test dR * ket_R ≈ ket_0
    @test dL' * ket_0 ≈ ket_L
    @test dR' * ket_0 ≈ ket_R
    @test dL' * ket_R ≈ -ket_LR
    @test dR' * ket_L ≈ ket_LR
    @test dR' * dL' * ket_0 ≈ ket_LR
end

@testset "Test DQD basis and operators in the g-e basis (no-blockade)" begin
    dqd = Dqd(0.5, 0.1, 1., 0., 0., 0.)
    ket_0, ket_g, ket_e, ket_ge = build_dqd_basis_ge(dqd)
    dg, de = build_dqd_fermi_ops_ge(dqd)

    @test dg * ket_g ≈ ket_0
    @test de * ket_e ≈ ket_0
    @test dg' * ket_0 ≈ ket_g
    @test de' * ket_0 ≈ ket_e
    @test de' * ket_g ≈ ket_ge
    @test dg' * ket_e ≈ -ket_ge
    @test de' * dg' * ket_0 ≈ ket_ge
end

@testset "Test DQD basis and operators in the left-right basis (blockade)" begin
    dqd = Dqd(17.57, 0.0, 6.78, 0., 0.)
    ket_0, ket_L, ket_R = build_dqd_basis_LR(dqd)
    dL, dR = build_dqd_fermi_ops_LR(dqd)

    @test dL * ket_L ≈ ket_0
    @test dR * ket_R ≈ ket_0
    @test dL' * ket_0 ≈ ket_L
    @test dR' * ket_0 ≈ ket_R
end

@testset "Test DQD basis and operators in the ground-excited basis (blockade)" begin
    dqd = Dqd(17.57, 0.0, 6.78, 0., 0.)
    ket_0, ket_g, ket_e = build_dqd_basis_ge(dqd)
    dg, de = build_dqd_fermi_ops_ge(dqd)

    @test dg * ket_g ≈ ket_0
    @test de * ket_e ≈ ket_0
    @test dg' * ket_0 ≈ ket_g
    @test de' * ket_0 ≈ ket_e
end

@testset "Test DQD Hamiltonian in the L-R basis (blockade)" begin
    dqd = Dqd(0.5, 0.1, 1., 0., 0.)
    ket_0, ket_L, ket_R = build_dqd_basis_LR(dqd)
    H_dqd_LR = build_H_dqd_LR(dqd)
    ϵg, ϵe = get_eigen_energies(dqd)
    evals = eigenenergies(H_dqd_LR)

    @test ϵg ≈ minimum(evals)
    @test ϵe ≈ maximum(evals)
end

@testset "Test DQD Hamiltonian in the g-e basis (blockade)" begin
    dqd = Dqd(0.5, 0.1, 1., 0., 0.)
    ket_0, ket_g, ket_e = build_dqd_basis_ge(dqd)
    H_dqd_ge = build_H_dqd_ge(dqd)
    ϵg, ϵe = get_eigen_energies(dqd)
    θ = get_θ(dqd)
    
    @test rad2deg(θ) ≈ rad2deg(acos(- dqd.Δϵ / get_Ω(dqd)))
    @test rad2deg(θ) ≈ rad2deg(atan(2 * dqd.tc, dqd.Δϵ))
    @test H_dqd_ge * ket_g ≈ ϵg * ket_g
    @test H_dqd_ge * ket_e ≈ ϵe * ket_e
    @test H_dqd_ge ≈ ϵg * ket_g * ket_g' + ϵe * ket_e * ket_e'
    @test H_dqd_ge.data ≈ build_H_dqd_LR(dqd).data
end

@testset "Test DQD Hamiltonian in the g-e basis (no-blockade)" begin
    dqd = Dqd(1., 0.3, 1.1, 0., 0., 0.)
    ket_0, ket_g, ket_e, ket_ge = build_dqd_basis_ge(dqd)
    H_dqd_ge = build_H_dqd_ge(dqd)
    ϵg, ϵe = get_eigen_energies(dqd)

    @test H_dqd_ge * ket_g ≈ ϵg * ket_g
    @test H_dqd_ge * ket_e ≈ ϵe * ket_e
    @test H_dqd_ge ≈ ϵg * ket_g * ket_g' + ϵe * ket_e * ket_e'
    @test H_dqd_ge ≈ build_H_dqd_LR(dqd)
end