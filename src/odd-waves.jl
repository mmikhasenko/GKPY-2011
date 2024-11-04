# P-wave
function cotδ1_less_1050(s; pars)
    s < 4mπ^2 && return 0.0
    @unpack e0, B0, B1 = pars
    return sqrt(s) / (2 * kπ(s)^3) *
           (mρ^2 - s) *
           (2mπ^3 / (mρ^2 * sqrt(s)) + B0 + B1 * conformal_w(s; s0 = e0^2))
end
function δ1_more_1050(s; pars)
    s < 4mπ^2 && return 0.0
    @unpack λ1, λ2, ϵ1, ϵ2 = pars
    λ0 = acot(cotδ1_less_1050(4mK^2; pars))
    return λ0 + λ1 * (sqrt(s) / (2mK) - 1) + λ2 * (sqrt(s) / (2mK) - 1)^2
end
function _δ1(s; pars = p_wave_pars)
    s < 4mπ^2 && return 0.0
    @unpack e0 = pars
    if s < e0^2
        v = acot(cotδ1_less_1050(s; pars))
        return (v < 0) ? v + π : v
    end
    return s < 1.4^2 ? δ1_more_1050(s; pars) : δ1_more_1050(1.4^2; pars)
end

"""
    δ1(s; pars=p_wave_pars)

Computes the phase shift of the P-wave.
The default value of parameters are taken from the paper (see p_wave_pars).

# Example
```julia
julia> δ1(0.5^2)
0.7547425109417016
```
"""
function δ1(model::GKPY11, s; pars = model.P)
    v = _δ1(s; pars)
    s < 0.8^2 ? v : (v > 0 ? v : v + π)
end

"""
    cotδ1(s; pars=p_wave_pars)

Computes the cotangence of the P-wave phase shift.
The default value of parameters are taken from the paper (see GKPY11.p_wave_pars).
"""
cotδ1(model::GKPY11, s; pars = model.P) = cot(δ1(model, s; pars))

p_wave_phase_shift(model::GKPY11, s::Real; pars = model.P) = δ1(model, s; pars)
p_wave_elasticity(model::GKPY11, s::Real; pars = model.P) = one(s)

"""
    p_wave_amplitude(model::GKPY11, s)

Computes the P-wave scattering amplitude for a given model `model` at energy squared `s` (in GeV²).
The amplitude is computed from the P-wave phase shift and elasticity, and returns a complex value.

## See also:
- [`p_wave_phase_shift`](@ref p_wave_phase_shift): Computes the P-wave phase shift
- [`p_wave_elasticity`](@ref p_wave_elasticity): Computes the P-wave elasticity.
"""
function p_wave_amplitude(model::GKPY11, s::Real; pars = model.P)
    _δ = p_wave_phase_shift(model, s; pars)
    _η = p_wave_elasticity(model, s; pars)
    _t = amplitude_from_phase_and_elasticity(_δ, _η, mπ)
    return _t
end

# F-wave
"""
    cotδ3(model::GKPY11, s; pars = model.F)

Computes the cotangence of the F-wave phase shift.
The default value of parameters are taken from the paper (see GKPY11.f_wave_pars).
"""
function cotδ3(model::GKPY11, s; pars = model.F)
    @unpack B0, B1, λ = pars
    w = conformal_w(s; s0 = 1.45^2)
    return sqrt(s) / (2 * kπ(s)^7) * mπ^6 * (2λ * mπ / sqrt(s) + B0 + B1 * w)
end

δ3(model::GKPY11, s; pars = model.F) = acot(cotδ3(model, s; pars))

"""
    f_wave_phase_shift(model::GKPY11, s::Real)

Computes the phase shift of the F-wave.
The default value of parameters are taken from the paper (see GKPY11.f_wave_pars).

# Example
```julia
julia> f_wave_phase_shift(GKPY11.f_default, 0.5^2)
0.04172757888407087
```
"""
f_wave_phase_shift(model::GKPY11, s::Real; pars = model.F) = δ3(model, s; pars)
f_wave_elasticity(model::GKPY11, s::Real; pars = model.F) = one(s)

"""
    f_wave_amplitude(model::GKPY11, s; pars = model.F)

Computes the P-wave scattering amplitude for a given model `model` at energy squared `s` (in GeV²).
The amplitude is computed from the P-wave phase shift and elasticity, and returns a complex value.

## See also:
- [`f_wave_phase_shift`](@ref f_wave_phase_shift): Computes the P-wave phase shift in radians.
- [`f_wave_elasticity`](@ref f_wave_elasticity): Computes the P-wave elasticity.
"""
function f_wave_amplitude(model::GKPY11, s::Real; pars = model.F)
    _δ = f_wave_phase_shift(model, s; pars)
    _η = f_wave_elasticity(model, s; pars)
    _t = amplitude_from_phase_and_elasticity(_δ, _η, mπ)
    return _t
end
