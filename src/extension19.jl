function conformal_f_KK(s, coeff)
    w1 = conformal_w1(s; s_left = (2mK)^2, s_right = 1.5^2)
    _sum = sum(enumerate(coeff)) do (i, Ki)
        xi = CheybyshevT{i - 1}()
        Ki * xi(w1)
    end
    return mK^2 * _sum
end

"""
    d_f0(s, M, G, coeff)

Calculate the denominator of the t_f0 amplitude.

## Arguments
  * `s::Complex` : Mandelstam variable
  * `Msq::Real` : mass parameter
  * `G::Real` : width parameter
  * `coeff::Tuple` : coefficients of the K polynomial
"""
d_f0(s, Msq, G, coeff) =
    Msq - s - J_bar(s, mπ) * s * G - J_bar(s, mK) * conformal_f_KK(s, coeff)

function t_f0(s::Complex; pars)
    @unpack sp, coeff_K = pars
    @unpack Msq, G = get_MG(sp, coeff_K)
    return s * G / d_f0(s, Msq, G, coeff_K)
end
t_f0(s::Real; pars) = t_f0(s + iϵ; pars)


function get_MG(s_p, coeff_K)
    #
    pol1(Msq, G) = d_f0(s_p, Msq, G, coeff_K) + 2im * σ(s_p, mπ) * s_p * G
    #
    C = pol1(0, 0) # constant term
    #
    Δpol1(Msq, G) = pol1(Msq, G) - C
    A, B = Δpol1(1, 0), Δpol1(0, 1) # first-order terms
    #
    X = [reim(A) |> collect reim(B) |> collect] # matrix, first row - real, second - imag
    b = (reim(-C) |> collect) # vector of the left part
    # invert matrix
    Msq, G = X \ b
    #
    (; Msq, G)
end

struct GKPY11_pure <: GKPY11
end
t_conf(s; pars) = 1 / σ(s, mπ) / (cotδ0_interval1(GKPY11_pure(), s; pars) - 1im)

## Interval 1


"""
    t0_interval1(s; pars)

The scattering amplitude of the interval1 is computed from a combination of the conformal amplitude and the f0 amplitude.
The two parametrizations are combined on the level of S-matrix elements, S = S_conf * S_f0.
"""
function t0_interval1(model::PRR19, s; pars = model.S)
    _tconf = t_conf(s; pars)
    _tf0 = t_f0(s; pars)
    return _tconf + _tf0 + 2im * σ(s, mπ) * _tconf * _tf0
end

function δ0_interval1_deg(model::PRR19, s; pars = model.S)
    _t0 = t0_interval1(model, s; pars)
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


"""
    η0_interval1(s; pars)

Inelasticity for the first interval is computed from the scattering amplitude,
by taking absolute value of `η = |1 + 2iσ t|`.
"""
η0_interval1(model::PRR19, s; pars = model.S) = abs(1 + t0_interval1(model, s; pars) * 2im * σ(s, mπ))



## Interval 2

"""
    δ0_interval2_deg(s; pars)

Calculate the phase shift for the second interval, 1.4 < sqrt(s) < 2.0 GeV.
It extends the phase shift from the first interval, by polynomial dependence.
"""
function δ0_interval2_deg(model::PRR19, s::Real; pars = model.S)
    @unpack coeff_d = pars
    #
    m_boundary1 = 1.4
    sm = m_boundary1^2
    w2(sx) = conformal_w1(sx; s_left = sm, s_right = 2.0^2)
    δ0_boundary1 = δ0_interval1_deg(model, sm; pars)
    ϵ = 1e-6
    dδ0_ds = (δ0_interval1_deg(model, sm + ϵ; pars) - δ0_boundary1) / ϵ
    dw2_ds = (w2(sm + ϵ) - w2(sm - ϵ)) / (2ϵ)
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

δ0_interval2(model::PRR19, s; pars = model.S) = δ0_interval2_deg(model, s; pars) / 180 * π

"""
    t0_interval2(s; pars)

Calculate the scattering amplitude in the real axis for the second interval,
1.4 < sqrt(s) < 2.0 GeV, using phenomenological parametrization of the phase shift and inelasticity.
"""
function t0_interval2(model::PRR19, s::Real; pars = model.S)
    _δ0 = δ0_interval2(model, s; pars)
    _eta = η0_interval2(model, s; pars)
    return (exp(2im * _δ0) * _eta - 1) / (2im * σ(s, mπ))
end


"""
    η0_interval2(s; pars)

Inelasticity for the second interval is computed using phenomenological parametrization.
It is a exponent of a square of a polynomial of a break-up momentum.
"""
function η0_interval2(model::PRR19, s; pars = model.S)
    @unpack coeff_ϵ_from2, m_boundary1 = pars
    sm = m_boundary1^2
    #
    _ηm = η0_interval1(model, sm; pars)
    ϵ0 = sqrt(-log(_ηm))
    ϵ = 1e-6
    dη0_ds = (η0_interval1(model, sm + ϵ; pars) - _ηm) / ϵ
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

function δ0_deg(model::PRR19, s::Real; pars = model.S)
    @unpack m_boundary1, m_boundary2 = pars
    s < m_boundary1^2 && return δ0_interval1_deg(model, s; pars)
    s < m_boundary2^2 && return δ0_interval2_deg(model, s; pars)
    δ0_interval2_deg(model, m_boundary2^2; pars)
end

s_wave_phase_shift(model::PRR19, s::Real; pars = model.S) = δ0_deg(model, s; pars) / 180 * π

function s_wave_elasticity(model::PRR19, s::Real; pars = model.S)
    @unpack m_boundary1, m_boundary2 = pars
    s < m_boundary1^2 && return η0_interval1(model, s; pars)
    s < m_boundary2^2 && return η0_interval2(model, s; pars)
    η0_interval2(model, m_boundary2^2; pars)
end

function s_wave_amplitude(model::PRR19, s::Real; pars = model.S)
    _δ = s_wave_phase_shift(model, s; pars)
    _η = s_wave_elasticity(model, s; pars)
    _t = amplitude_from_phase_and_elasticity(_δ, _η, σ(s, mπ))
    return _t
end
