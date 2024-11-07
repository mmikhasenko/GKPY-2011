@with_kw struct GKPY11{tS, tP, tD, tF} <: PiPiBase
    S::tS = NamedTuple()
    P::tP = NamedTuple()
    D::tD = NamedTuple()
    F::tF = NamedTuple()
end

GKPY11_default = GKPY11(;
    S = (m_boundary1 = 0.85,
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
    ),
    P = (B0 = 1.043, B1 = 0.19, λ1 = 1.39, λ2 = -1.70, ϵ1 = 0.00, ϵ2 = 0.07, e0 = 1.05),
    D = NamedTuple(),
    F = (B0 = 1.09e5, B1 = 1.41e5, λ = 0.051e5),
)
