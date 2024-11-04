const mπ = 0.13957
const mρ = 0.7736
const Γρ = 0.146
const mK = 0.496
const mη = 0.54751
#
const iϵ = 1e-6im

breakup(s, m) = sqrt(s / 4 - m^2)
σ(s, m) = 2 * breakup(s, m) / sqrt(s)
J_bar(s, m) = 2 / π + σ(s, m) / π * log((σ(s, m) - 1) / (σ(s, m) + 1))

conformal_w(s; s0::Float64) = (sqrt(s) - sqrt(s0 - s)) / (sqrt(s) + sqrt(s0 - s))

conformal_w1(s; s_left::Float64, s_right::Float64) =
    2 * (sqrt(s) - sqrt(s_left)) / (sqrt(s_right) - sqrt(s_left)) - 1

kπ(s) = breakup(s, mπ)
kK(s) = breakup(s, mK)
kη(s) = breakup(s, mη)

amplitude_from_phase_and_elasticity(δ, η, m) = (η .* exp(2im * δ) - 1) / (2σ(s, m))
