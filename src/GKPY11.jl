module GKPY11

using Parameters

import ..PiPiMadrid: iϵ, mπ, mK, mη, mρ, Γρ
import ..PiPiMadrid: breakup, σ, conformal_w

kπ(s) = breakup(s, mπ)
kK(s) = breakup(s, mK)
kη(s) = breakup(s, mη)

include("odd-waves.jl")

include("s-wave.jl")

end
