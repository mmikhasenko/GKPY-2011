using Test
using PiPiMadrid.Parameters
using PiPiMadrid

@testset "Test δ1 function" begin
    @test p_wave_phase_shift(GKPY11_default, 0.1) ≈ 0.005348662528080389
    @test p_wave_phase_shift(GKPY11_default, 0.5) ≈ 0.7547425109417016
    @test p_wave_phase_shift(GKPY11_default, 1.5) ≈ 2.9163793203966426
    @test p_wave_phase_shift(GKPY11_default, 2.0) ≈ 2.967957560740781
end

@testset "Test δ3 function" begin
    @test f_wave_phase_shift(GKPY11_default, 0.1) ≈ 4.4539237752958584e-7
    @test f_wave_phase_shift(GKPY11_default, 0.5) ≈ 0.0020545780665097556
    @test f_wave_phase_shift(GKPY11_default, 1.5) ≈ 0.04172757888407087
    @test f_wave_phase_shift(GKPY11_default, 2.0) ≈ 0.07382739896763756
end

@testset "Test δ0 function" begin
    @test s_wave_phase_shift(GKPY11_default, 0.5^2) ≈ 39.86680197881492 / 180 * π
    @test s_wave_phase_shift(GKPY11_default, 1.1^2) ≈ 249.68337945026678 / 180 * π
    @test s_wave_phase_shift(GKPY11_default, 1.4^2) ≈ 310.8928743137025 / 180 * π
    @test s_wave_phase_shift(GKPY11_default, 2.0^2) ≈ 359.88122836077986 / 180 * π
end


# Internal functions

import PiPiMadrid: cotδ1

@testset "Test cotδ1 function" begin
    @test PiPiMadrid.cotδ1(GKPY11_default, 0.1) ≈ 186.9608446351338
    @test PiPiMadrid.cotδ1(GKPY11_default, 0.5) ≈ 1.0632707316721013
    @test PiPiMadrid.cotδ1(GKPY11_default, 1.5) ≈ -4.36490825792332
    @test PiPiMadrid.cotδ1(GKPY11_default, 2.0) ≈ -5.701209421847014
end

import PiPiMadrid: cotδ3

@testset "Test cotδ3 function" begin
    @test PiPiMadrid.cotδ3(GKPY11_default, 0.1) ≈ 2.2452113023274797e6
    @test PiPiMadrid.cotδ3(GKPY11_default, 0.5) ≈ 486.7172531445134
    @test PiPiMadrid.cotδ3(GKPY11_default, 1.5) ≈ 23.951054971265382
    @test PiPiMadrid.cotδ3(GKPY11_default, 2.0) ≈ 13.520488669013462
end

import PiPiMadrid: CheybyshevT
@testset "Test CheybyshevT function" begin
    @test PiPiMadrid.CheybyshevT{0}()(0.5) ≈ 1.0
    @test PiPiMadrid.CheybyshevT{1}()(0.5) ≈ 0.5
    @test PiPiMadrid.CheybyshevT{2}()(0.5) ≈ -0.5
    @test PiPiMadrid.CheybyshevT{3}()(0.5) ≈ -1.0
    @test PiPiMadrid.CheybyshevT{4}()(0.5) ≈ -0.5
    @test PiPiMadrid.CheybyshevT{5}()(0.5) ≈ 0.5
end

import PiPiMadrid: iϵ, mπ, J_bar, σ

@testset "Test Chew-Mandelsta" begin
    @test isapprox(imag(J_bar(0.5^2 + iϵ, mπ)), σ(0.5^2, mπ), atol = 1e-5)
    @test isapprox(imag(J_bar(0.9^2 + iϵ, mπ)), σ(0.9^2, mπ), atol = 1e-5)
end

import PiPiMadrid: δ0_interval1_deg, δ0_interval2_deg
import PiPiMadrid: η0_interval1, η0_interval2
import PiPiMadrid: t0_interval1, t0_interval2

@testset "Test of matching point" begin
    @unpack m_boundary1 = PRR19_default.S
    sm = m_boundary1^2
    #
    @test PiPiMadrid.δ0_interval1_deg(PRR19_default, sm) ≈ PiPiMadrid.δ0_interval2_deg(PRR19_default, sm)
    @test PiPiMadrid.η0_interval1(PRR19_default, sm) ≈ PiPiMadrid.η0_interval2(PRR19_default, sm)
    @test PiPiMadrid.t0_interval1(PRR19_default, sm) ≈ PiPiMadrid.t0_interval2(PRR19_default, sm)
end

@testset "Elasticity interval 1" begin
    model = PiPiMadrid.PRR19_default
    s_test = 1.1^2
    @test PiPiMadrid.η0_interval1(model, s_test) ≈ 0.4872435393540931
    @test PiPiMadrid.t0_interval1(model, s_test) ≈ 0.18573254487071456 + 0.6870385014284162im
end

@testset "Close to threshold" begin
    m_thr = m_thr = 2mπ + 8.627191629750897e-6
    v = PiPiMadrid.t0_interval1(PRR19_default, m_thr^2)
    # Imaginary part should be almost zero
    @test imag(v) ≈ 0.00041310218334773683
    @test real(v) ≈ 0.22823338085075165
end
