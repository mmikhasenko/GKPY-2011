abstract type GKPY11 <: PiPiBase end

@with_kw struct GKPY11_default{tS, tP, tD, tF} <: GKPY11
    S::tS = (m_boundary1 = 0.85,
        m_boundary3 = 1.42,
        # interval 1,
        b_coeffs = (7.14, -25.3, -33.2, -26.2),
        z0 = 0.13957,
        # interval 1,
        d0 = 226.5,
        c = -81.0,
        # interval 2,
        b = 93.3,
        c_coeff = 48.7,
        d = -88.3,
    )
    P::tP = (B0 = 1.043, B1 = 0.19, λ1 = 1.39, λ2 = -1.70, ϵ1 = 0.00, ϵ2 = 0.07, e0 = 1.05)
    D::tD = NamedTuple()
    F::tF = (B0 = 1.09e5, B1 = 1.41e5, λ = 0.051e5)
end
