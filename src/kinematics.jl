# masses from the paper
const mπ = 0.13957
const mρ = 0.7736
const Γρ = 0.146
const mK = 0.496
const mη = 0.54751

# maps
breakup(s, m) = sqrt(s / 4 - m^2)

kπ(s) = breakup(s, mπ)
kK(s) = breakup(s, mK)
kη(s) = breakup(s, mη)

conformal_w(s; s0::Float64) = (sqrt(s) - sqrt(s0 - s)) / (sqrt(s) + sqrt(s0 - s))

const iϵ = 1e-6im
