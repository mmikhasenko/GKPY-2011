@with_kw struct PRR19{tS} <: PiPiBase
    S::tS = NamedTuple()
end

const PRR19_default = PRR19(
    # default parameters from the paper
    S = (
        b_coeffs = (12.2, -0.9, 15.9, -5.7, -22.5, 6.9),
        z0 = 0.137, # Adler zero at z0^2/2
        coeff_K = (5.25, -4.4, 0.175, -0.28),
        sp = (0.996 + 0.025im)^2,
        #
        m_boundary1 = 1.4,
        coeff_d = (-5.4, 0.0, 0.0),
        coeff_ϵ_from2 = (10.3, 0.0, 0.0),
        #
        m_boundary2 = 2.0,
    ),
)
