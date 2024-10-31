
@with_kw struct SWavePars
    m_boundary1::Float64 = 0.85
    m_boundary3::Float64 = 1.42
    # interval 1
    b_coeffs::NTuple{4, Float64} = (7.14, -25.3, -33.2, -26.2)
    z0::Float64 = 0.13957 # Adler zero at z0^2/2
    # interval 1
    d0::Float64 = 226.5
    c::Float64 = -81.0
    # interval 2
    b::Float64 = 93.3
    c_coeff::Float64 = 48.7
    d::Float64 = -88.3
end

const s_wave_pars = SWavePars();

"""
    δ0(s; pars=s_wave_pars)

Computes the phase shift of the F-wave.
The default value of parameters are taken from the paper (see GKPY11.s_wave_pars).

# Example
```julia
julia> δ0(0.5^2)
39.86680197881492
```
"""
function δ0(s::Real; pars = s_wave_pars)
    @unpack m_boundary1, m_boundary3 = pars
    if s <= (2mπ)^2
        return 0.0
    elseif s <= (m_boundary1)^2
        return δ0_interval1(s; pars)
    elseif s <= (2 * mK)^2
        return δ0_interval2(s; pars)
    elseif s <= (m_boundary3)^2
        return δ0_interval3(s; pars)
    end
    return δ0_inverval4(s; pars)
end


function cotδ0_interval1(s::Complex; pars = s_wave_pars)
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

cotδ0_interval1(s::Real; pars = s_wave_pars) = cotδ0_interval1(s + iϵ; pars) |> real

function δ0_interval1(s; pars = s_wave_pars)
    _cotδ0_interval1 = cotδ0_interval1(s; pars)
    _δ0_interval1 = atan(1 / _cotδ0_interval1)
    return _δ0_interval1 < 0 ? _δ0_interval1 + π : _δ0_interval1
end

function δ0_interval2(s; pars = s_wave_pars)
    @unpack m_boundary1, d0, c = pars
    #
    ϵ = 1e-6
    δ1_boundary1 = δ0_interval1(m_boundary1^2; pars)
    dδ1_dm_boundary1 =
        (δ0_interval1(m_boundary1^2 + ϵ; pars) - δ0_interval1(m_boundary1^2 - ϵ; pars)) /
        (2ϵ)
    #
    δ1_boundary1_deg = δ1_boundary1 / π * 180
    dδ1_dm_boundary1_deg = dδ1_dm_boundary1 / π * 180
    #
    a_k2 = s > (2 * mK)^2 ? kK(s) : sqrt(mK^2 - s / 4)
    k2_m = sqrt(mK^2 - m_boundary1^2 / 4)
    _δ0 = d0 * (1 - a_k2 / k2_m)^2 +
          δ1_boundary1_deg * a_k2 / k2_m * (2 - a_k2 / k2_m) +
          a_k2 * (k2_m - a_k2) * (8 * dδ1_dm_boundary1_deg + c * (k2_m - a_k2) / (mK^3))
    return _δ0 / 180.0 * π
end

function δ0_interval3(s; pars = s_wave_pars)
    @unpack d0, b, d, c_coeff = pars
    k2 = kK(s)
    _δ0 = d0 +
          b * (k2 / mK)^2 +
          c_coeff * (k2 / mK)^4 +
          d * (s > (2 * mη)^2 ? (kη(s) / mη)^2 : 0)
    return _δ0 / 180.0 * π
end

function δ0_inverval4(s; pars = s_wave_pars)
    @unpack m_boundary3 = s_wave_pars
    #
    ϵ = 1e-6
    δ3_boundary3 = δ0_interval3(m_boundary3^2; pars)
    dδ3_dm_boundary3 =
        (δ0_interval3(m_boundary3^2 + ϵ; pars) - δ0_interval3(m_boundary3^2 - ϵ; pars)) /
        (2e-6)
    #
    a = 2π
    alpha = dδ3_dm_boundary3 / (a - δ3_boundary3)
    p = dδ3_dm_boundary3 / (a * alpha) * exp(alpha * (m_boundary3)^2)
    return a * (1 - p * exp(-alpha * s))
end
