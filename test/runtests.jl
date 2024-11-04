using Test
using PiPiMadrid.Parameters

using PiPiMadrid#
import PiPiMadrid: iϵ, mπ
import PiPiMadrid: J_bar, σ

@testset "Test δ1 function" begin
    @test p_wave_phase_shift(GKPY11_default(), 0.1) ≈ 0.005348662528080389
    @test p_wave_phase_shift(GKPY11_default(), 0.5) ≈ 0.7547425109417016
    @test p_wave_phase_shift(GKPY11_default(), 1.5) ≈ 2.9163793203966426
    @test p_wave_phase_shift(GKPY11_default(), 2.0) ≈ 2.967957560740781
end

@testset "Test cotδ1 function" begin
    @test PiPiMadrid.cotδ1(GKPY11_default(), 0.1) ≈ 186.9608446351338
    @test PiPiMadrid.cotδ1(GKPY11_default(), 0.5) ≈ 1.0632707316721013
    @test PiPiMadrid.cotδ1(GKPY11_default(), 1.5) ≈ -4.36490825792332
    @test PiPiMadrid.cotδ1(GKPY11_default(), 2.0) ≈ -5.701209421847014
end

@testset "Test δ3 function" begin
    @test f_wave_phase_shift(GKPY11_default(), 0.1) ≈ 4.4539237752958584e-7
    @test f_wave_phase_shift(GKPY11_default(), 0.5) ≈ 0.0020545780665097556
    @test f_wave_phase_shift(GKPY11_default(), 1.5) ≈ 0.04172757888407087
    @test f_wave_phase_shift(GKPY11_default(), 2.0) ≈ 0.07382739896763756
end

@testset "Test cotδ3 function" begin
    @test PiPiMadrid.cotδ3(GKPY11_default(), 0.1) ≈ 2.2452113023274797e6
    @test PiPiMadrid.cotδ3(GKPY11_default(), 0.5) ≈ 486.7172531445134
    @test PiPiMadrid.cotδ3(GKPY11_default(), 1.5) ≈ 23.951054971265382
    @test PiPiMadrid.cotδ3(GKPY11_default(), 2.0) ≈ 13.520488669013462
end

@testset "Test δ0 function" begin
    @test s_wave_phase_shift(GKPY11_default(), 0.5^2) ≈ 39.86680197881492 / 180 * π
    @test s_wave_phase_shift(GKPY11_default(), 1.1^2) ≈ 249.68337945026678 / 180 * π
    @test s_wave_phase_shift(GKPY11_default(), 1.4^2) ≈ 310.8928743137025 / 180 * π
    @test s_wave_phase_shift(GKPY11_default(), 2.0^2) ≈ 359.88122836077986 / 180 * π
end

@testset "Test CheybyshevT function" begin
    @test PiPiMadrid.CheybyshevT{0}()(0.5) ≈ 1.0
    @test PiPiMadrid.CheybyshevT{1}()(0.5) ≈ 0.5
    @test PiPiMadrid.CheybyshevT{2}()(0.5) ≈ -0.5
    @test PiPiMadrid.CheybyshevT{3}()(0.5) ≈ -1.0
    @test PiPiMadrid.CheybyshevT{4}()(0.5) ≈ -0.5
    @test PiPiMadrid.CheybyshevT{5}()(0.5) ≈ 0.5
end

@testset "Test Chew-Mandelsta" begin
    @test isapprox(imag(J_bar(0.5^2 + iϵ, mπ)), σ(0.5^2, mπ), atol = 1e-5)
    @test isapprox(imag(J_bar(0.9^2 + iϵ, mπ)), σ(0.9^2, mπ), atol = 1e-5)
end

# @testset "Test conformal amplitude" begin
#     @test PRR19.t_conf(0.25; pars = GKPY11.s_wave_pars) ≈ 0.5930131073223124 + 0.49525238695313145im
# end

# @testset "Test f0 amplitude" begin
#     @unpack sp = PRR19.s_wave_f0_pars
#     @test PRR19.t_f0(sp; pars = PRR19.s_wave_f0_pars) * 2im * σ(sp, mπ) ≈ -1.0
#     @test PRR19.t_f0(sp'; pars = PRR19.s_wave_f0_pars) * 2im * σ(sp', mπ) ≈ 1.0
# end

# @testset "Test of matching point" begin
#     pars = PRR19.s_wave_f0_pars
#     #
#     m_boundary1 = 1.4
#     sm = m_boundary1^2
#     @test PRR19.δ0_interval1_deg(sm; pars) ≈ PRR19.δ0_interval2_deg(sm; pars)
#     @test PRR19.η0_interval1(sm; pars) ≈ PRR19.η0_interval2(sm; pars)
#     @test PRR19.t0_interval1(sm; pars) ≈ PRR19.t0_interval2(sm; pars)
# end
