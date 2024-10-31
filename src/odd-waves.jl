# fit parameters (constrained)
const p_wave_pars = (B0 = 1.043, B1 = 0.19, λ1 = 1.39, λ2 = -1.70, ϵ1 = 0.00, ϵ2 = 0.07, e0 = 1.05)
const f_wave_pars = (B0 = 1.09e5, B1 = 1.41e5, λ = 0.051e5)


# P-wave
function cotδ1_less_1050(s; pars)
    s < 4mπ^2 && return 0.0
    @unpack e0, B0, B1 = pars
    return sqrt(s) / (2 * kπ(s)^3) * (mρ^2 - s) * (2mπ^3 / (mρ^2 * sqrt(s)) + B0 + B1 * conformal_w(s; s0 = e0^2))
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
The default value of parameters are taken from the paper (see GKPY11.p_wave_pars).

# Example
```julia
julia> δ1(0.5^2)
0.7547425109417016
```
"""
function δ1(s; pars = p_wave_pars)
    v = _δ1(s; pars = p_wave_pars)
    s < 0.8^2 ? v : (v > 0 ? v : v + π)
end

"""
    cotδ1(s; pars=p_wave_pars)

Computes the cotangence of the P-wave phase shift.
The default value of parameters are taken from the paper (see GKPY11.p_wave_pars).
"""
cotδ1(s; pars = p_wave_pars) = cot(δ1(s; pars))

# D-wave
"""
    cotδ3(s; pars=f_wave_pars)

Computes the cotangence of the F-wave phase shift.
The default value of parameters are taken from the paper (see GKPY11.f_wave_pars).
"""
function cotδ3(s; pars = f_wave_pars)
    @unpack B0, B1, λ = pars
    w = conformal_w(s; s0 = 1.45^2)
    return sqrt(s) / (2 * kπ(s)^7) * mπ^6 * (2λ * mπ / sqrt(s) + B0 + B1 * w)
end

"""
    δ3(s; pars=f_wave_pars)

Computes the phase shift of the F-wave.
The default value of parameters are taken from the paper (see GKPY11.f_wave_pars).

# Example
```julia
julia> δ3(0.5^2)
0.04172757888407087
```
"""
δ3(s; pars = f_wave_pars) = acot(cotδ3(s; pars))
