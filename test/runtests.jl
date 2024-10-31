using Test
using GKPY11

@testset "Test δ1 function" begin
    @test δ1(0.1) ≈ 0.005348662528080389
    @test δ1(0.5) ≈ 0.7547425109417016
    @test δ1(1.5) ≈ 2.9163793203966426
    @test δ1(2.0) ≈ 2.967957560740781
end

@testset "Test cotδ1 function" begin
    @test cotδ1(0.1) ≈ 186.9608446351338
    @test cotδ1(0.5) ≈ 1.0632707316721013
    @test cotδ1(1.5) ≈ -4.36490825792332
    @test cotδ1(2.0) ≈ -5.701209421847014
end

@testset "Test δ3 function" begin
    @test δ3(0.1) ≈ 4.4539237752958584e-7
    @test δ3(0.5) ≈ 0.0020545780665097556
    @test δ3(1.5) ≈ 0.04172757888407087
    @test δ3(2.0) ≈ 0.07382739896763756
end

@testset "Test cotδ3 function" begin
    @test cotδ3(0.1) ≈ 2.2452113023274797e6
    @test cotδ3(0.5) ≈ 486.7172531445134
    @test cotδ3(1.5) ≈ 23.951054971265382
    @test cotδ3(2.0) ≈ 13.520488669013462
end

@testset "Test δ0 function" begin
    @test δ0(0.5^2) ≈ 39.86680197881492 / 180 * π
    @test δ0(1.1^2) ≈ 249.68337945026678 / 180 * π
    @test δ0(1.4^2) ≈ 310.8928743137025 / 180 * π
    @test δ0(2.0^2) ≈ 359.88122836077986 / 180 * π
end



import GKPY11: mπ, mK, breakup, cotδ0_interval1, s_wave_pars, iϵ
using Parameters
using Plots


struct CheybyshevT{N} end
(ch::CheybyshevT{0})(x::T) where T = 1.0
(ch::CheybyshevT{1})(x::T) where T = x
(ch::CheybyshevT{2})(x::T) where T = 2x^2 - 1
(ch::CheybyshevT{3})(x::T) where T = 4x^3 - 3x
(ch::CheybyshevT{4})(x::T) where T = 8x^4 - 8x^2 + 1
(ch::CheybyshevT{5})(x::T) where T = 16x^5 - 20x^3 + 5x


# default parameters from the paper
s_wave_f0_pars = (
    b_coeffs = (12.2, -0.9, 15.9, -5.7, -22.5, 6.9),
    z0 = 0.137, # Adler zero at z0^2/2
    coeff_K = (5.25, -4.4, 0.175, -0.28),
    sp = (0.996 + 0.025im)^2,
    #
    m_boundary1 = 1.4,
    coeff_d = (-5.4, 0.0, 0.0),
    coeff_ϵ_from2 = (10.3, 0.0, 0.0),
)

σ(s, m) = 2 * breakup(s, m) / sqrt(s)
J_bar(s, m) = 2 / π + σ(s, m) / π * log((σ(s, m) - 1) / (σ(s, m) + 1))

@test isapprox(imag(J_bar(0.5^2 + iϵ, mπ)), σ(0.5^2, mπ), atol = 1e-5)

conformal_w1(s; s_left::Float64, s_right::Float64) =
    2 * (sqrt(s) - sqrt(s_left)) / (sqrt(s_right) - sqrt(s_left)) - 1

function conformal_f_KK(s, coeff)
    w1 = conformal_w1(s; s_left = (2mK)^2, s_right = 1.5^2)
    return mK^2 * sum(enumerate(coeff)) do (i, Ki)
        xi = CheybyshevT{i - 1}()
        Ki * xi(w1)
    end
end


d_f0(s, M, G, coeff) =
    M - s - J_bar(s, mπ) * s * G - J_bar(s, mK) * conformal_f_KK(s, coeff)

function t_f0(s::Complex; pars)
    @unpack sp, coeff_K = pars
    @unpack M, G = get_MG(sp, coeff_K)
    return s * G / d_f0(s, M, G, coeff_K)
end
t_f0(s::Real; pars) = t_f0(s + iϵ; pars)


function get_MG(s_p, coeff_K)
    #
    pol1(M, G) = d_f0(s_p, M, G, coeff_K) + 2im * σ(s_p, mπ) * s_p * G
    #
    C = pol1(0, 0) # constant term
    #
    Δpol1(M, G) = pol1(M, G) - C
    A, B = Δpol1(1, 0), Δpol1(0, 1) # first-order terms
    #
    X = [reim(A) |> collect reim(B) |> collect] # matrix, first row - real, second - imag
    b = (reim(-C) |> collect) # vector of the left part
    # invert matrix
    M, G = X \ b
    #
    (; M, G)
end

let
    @unpack sp = s_wave_f0_pars
    @test t_f0(sp; pars = s_wave_f0_pars) * 2im * σ(sp, mπ) ≈ -1.0
    @test t_f0(sp'; pars = s_wave_f0_pars) * 2im * σ(sp', mπ) ≈ 1.0
end


t_conf(s; pars) = 1 / σ(s, mπ) / (cotδ0_interval1(s; pars) - 1im)
@test t_conf(0.25; pars = s_wave_pars) isa Complex

function t0_interval1(s; pars)
    _tconf = t_conf(s; pars)
    _tf0 = t_f0(s; pars)
    return _tconf + _tf0 + 2im * σ(s, mπ) * _tconf * _tf0
end

let
    plot()
    plot!(e -> real(t0_interval1(e^2; pars = s_wave_f0_pars)), 0.3, 1.5)
    plot!(e -> imag(t0_interval1(e^2; pars = s_wave_f0_pars)), 0.3, 1.5)
end

function δ0_interval1_deg(s; pars)
    _t0 = t0_interval1(s; pars)
    _δ0 = angle(_t0 * 2im * σ(s, mπ) + 1) / 2
    _δ0_def = _δ0 / π * 180
    shift = 0
    # ponts =
    sqrt_s = sqrt(s)
    if sqrt_s < 1.2 && _δ0_def < 0.0
        shift = 180
    end
    if 0.9 < sqrt_s < 1.5
        shift = 180
    end
    if 1.2 < sqrt_s < 1.8 && _δ0_def < 0.0
        shift = 360
    end
    if 1.5 < sqrt_s
        shift = 180 * 2
    end
    return _δ0_def + shift
end

let
    plot()
    plot!(e -> δ0_interval1_deg(e^2; pars = s_wave_f0_pars), 0.3, 1.4)
    plot!(e -> δ0_interval2_deg(e^2; pars = s_wave_f0_pars), 1.4, 2.0)
    vline!([0.9, 1.2, 1.5, 1.8])
end

δ0_interval1(s; pars) = δ0_interval1_deg(s; pars) / 180 * π

# Define δ0(s) using the given formula
function δ0_interval2_deg(s; pars)
    @unpack coeff_d = pars
    #
    m_boundary1 = 1.4
    w2(sx) = conformal_w1(sx; s_left = m_boundary1^2, s_right = 2.0^2)
    δ0_boundary1 = δ0_interval1_deg(m_boundary1^2; pars)
    ϵ = 1e-6
    dδ0_ds = (δ0_interval1_deg(m_boundary1^2 + ϵ; pars) -
              δ0_boundary1) / (ϵ)
    dw2_ds = (w2(m_boundary1^2 + ϵ) -
              w2(m_boundary1^2 - ϵ)) / (2ϵ)
    #
    _w2 = w2(s)
    d0, d1, d2 = coeff_d
    Δ0 = dδ0_ds / dw2_ds + 4 * d0 - 9 * d1 + 16 * d2
    _sum = sum(enumerate((Δ0, coeff_d...))) do (i, di)
        T = CheybyshevT{i}()
        pm1 = 2 * mod(i, 2) - 1
        return di * (T(_w2) + pm1)
    end
    _δ0 = δ0_boundary1 + _sum
    return _δ0
end

δ0_interval2(s; pars) = δ0_interval2_deg(s; pars) / 180 * π

function t0_interval2(s; pars)
    _δ0 = δ0_interval2(s; pars)
    _eta = η0_interval2(s; pars)
    return (exp(2im * _δ0) * _eta - 1) / (2im * σ(s, mπ))
end

η0_interval1(s; pars) = abs(1 + t0_interval1(s; pars) * 2im * σ(s, mπ))

function η0_interval2(s; pars)
    @unpack coeff_ϵ_from2, m_boundary1 = pars
    sm = m_boundary1^2
    #
    _ηm = η0_interval1(sm; pars)
    ϵ0 = sqrt(-log(_ηm))
    ϵ = 1e-6
    dη0_ds = (η0_interval1(sm + ϵ; pars) - _ηm) / ϵ
    #
    qm = breakup(sm, mπ)
    ϵ1 = -4qm^2 / ϵ0 * dη0_ds / _ηm
    coeff_ϵ = (ϵ0, ϵ1, coeff_ϵ_from2...)
    #
    Q = breakup(s, mπ) / qm - 1
    _sum = sum(enumerate(coeff_ϵ)) do (i, ϵi)
        ϵi * Q^(i - 1)
    end
    return exp(-_sum^2) |> real
end

let
    plot()
    plot!(e -> η0_interval1(e^2 + 1e-6im; pars = s_wave_f0_pars), 0.3, 1.4)
    plot!(e -> η0_interval2(e^2 + 1e-6im; pars = s_wave_f0_pars), 1.4, 2.0)
end

let
    plot()
    ev = range(0.3, 1.4, 1000)
    calv = map(e -> t0_interval1(e^2; pars = s_wave_f0_pars), ev)
    plot!(calv)
    ev = range(1.4, 2.0, 1000)
    calv = map(e -> t0_interval2(e^2; pars = s_wave_f0_pars), ev)
    plot!(calv)
end

let
    m_boundary1 = 1.4
    sm = m_boundary1^2
    t0_interval1(sm; pars = s_wave_f0_pars) - t0_interval2(sm; pars = s_wave_f0_pars)
end

let
    sm = 1.5
    pars = s_wave_f0_pars
    _t = t0_interval1(sm; pars)
    _δ = δ0_interval1(sm; pars)
    _η = η0_interval1(sm; pars)
    _t - (_η * exp(2im * _δ) - 1) / (2im * σ(sm, mπ))
    # (exp(2im * _δ0) * _eta - 1) / (2im * σ(s, mπ))
end

let
    m_boundary1 = 1.4
    sm = m_boundary1^2
    δ0_interval1_deg(sm; pars = s_wave_f0_pars) ≈ δ0_interval2_deg(sm; pars = s_wave_f0_pars)
end

let
    m_boundary1 = 1.4
    sm = m_boundary1^2
    pars = s_wave_f0_pars
    η0_interval1(sm; pars) ≈ η0_interval2(sm; pars)
end
