"""
    s_wave_phase_shift(model::GKPY11, s::Real; pars = model.S)

Computes the phase shift of the F-wave.
The default value of parameters are taken from the paper (see `GKPY11_default``).

# Example
```julia
julia> s_wave_phase_shift(0.5^2)
39.86680197881492
```
"""
function s_wave_phase_shift(model::GKPY11, s::Real; pars = model.S)
    @unpack m_boundary1, m_boundary3 = pars
    if s <= (2mπ)^2
        return 0.0
    elseif s <= (m_boundary1)^2
        return δ0_interval1(model, s; pars)
    elseif s <= (2 * mK)^2
        return δ0_interval2(model, s; pars)
    elseif s <= (m_boundary3)^2
        return δ0_interval3(model, s; pars)
    end
    return δ0_inverval4(model, s; pars)
end


function cotδ0_interval1(model::GKPY11, s::Complex; pars = model.S)
    @unpack b_coeffs, z0 = pars
    #
    w = conformal_w(s; s0 = (2 * mK)^2)
    conf_expansion = sum(enumerate(b_coeffs)) do (i, b)
        b * w^(i - 1)
    end
    _cotδ0_interval1 =
        sqrt(s) / (2 * kπ(s)) * mπ^2 / (s - z0^2 / 2) *
        (z0^2 / (mπ * sqrt(s)) + conf_expansion)
    return _cotδ0_interval1
end

cotδ0_interval1(model::GKPY11, s::Real; pars = model.S) =
    cotδ0_interval1(model, s + iϵ; pars)

function δ0_interval1(model::GKPY11, s::Real; pars = model.S)
    _cotδ0_interval1 = cotδ0_interval1(model, s; pars)
    _δ0_interval1 = atan(1 / _cotδ0_interval1)
    return _δ0_interval1 < 0 ? _δ0_interval1 + π : _δ0_interval1
end

function δ0_interval2(model::GKPY11, s::Real; pars = model.S)
    @unpack m_boundary1, d0, c = pars
    #
    ϵ = 1e-6
    _δ0_boundary1 = δ0_interval1(model, m_boundary1^2; pars)
    dδ0_dm_boundary1 =
        (δ0_interval1(model, m_boundary1^2 + ϵ; pars) - δ0_interval1(model, m_boundary1^2 - ϵ; pars)) /
        (2ϵ)
    #
    _δ0_boundary1_deg = _δ0_boundary1 / π * 180
    dδ0_dm_boundary1_deg = dδ0_dm_boundary1 / π * 180
    #
    a_k2 = s > (2 * mK)^2 ? kK(s) : sqrt(mK^2 - s / 4)
    k2_m = sqrt(mK^2 - m_boundary1^2 / 4)
    _δ0 = d0 * (1 - a_k2 / k2_m)^2 +
          _δ0_boundary1_deg * a_k2 / k2_m * (2 - a_k2 / k2_m) +
          a_k2 * (k2_m - a_k2) * (8 * dδ0_dm_boundary1_deg + c * (k2_m - a_k2) / (mK^3))
    return _δ0 / 180.0 * π
end

function δ0_interval3(model::GKPY11, s::Real; pars = model.S)
    @unpack d0, b, d, c_coeff = pars
    k2 = kK(s)
    _δ0 = d0 +
          b * (k2 / mK)^2 +
          c_coeff * (k2 / mK)^4 +
          d * (s > (2 * mη)^2 ? (kη(s) / mη)^2 : 0)
    return _δ0 / 180.0 * π
end

function δ0_inverval4(model::GKPY11, s::Real; pars = model.S)
    @unpack m_boundary3 = pars
    #
    ϵ = 1e-6
    δ3_boundary3 = δ0_interval3(model, m_boundary3^2; pars)
    dδ3_dm_boundary3 =
        (δ0_interval3(model, m_boundary3^2 + ϵ; pars) - δ0_interval3(model, m_boundary3^2 - ϵ; pars)) /
        (2e-6)
    #
    a = 2π
    alpha = dδ3_dm_boundary3 / (a - δ3_boundary3)
    p = dδ3_dm_boundary3 / (a * alpha) * exp(alpha * (m_boundary3)^2)
    return a * (1 - p * exp(-alpha * s))
end


"""
    s_wave_amplitude(pars::GKPY11, s; pars = model.p)

Computes the S-wave scattering amplitude for a given model `model` at energy squared `s` (in GeV²).
The amplitude is computed from the S-wave phase shift and elasticity, and returns a complex value.
The parameters, `pars` provided as an key argument overrides the model papametes

## See also:
- [`s_wave_phase_shift`](@ref s_wave_phase_shift): Computes the S-wave phase shift
- [`s_wave_elasticity`](@ref s_wave_elasticity): Computes the S-wave elasticity.
"""
function s_wave_amplitude(model::GKPY11, s::Real; pars = model.S)
    _δ = s_wave_phase_shift(model, s; pars)
    _η = s_wave_elasticity(model, s; pars)
    _t = amplitude_from_phase_and_elasticity(_δ, _η, σ(s, mπ))
    return _t
end
