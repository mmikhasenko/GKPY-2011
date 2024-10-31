module PRR19
using Parameters

import GKPY11: mπ, mK, breakup, cotδ0_interval1, iϵ

export t0_interval1, t0_interval2, t0
export δ0_interval1, δ0_interval2, δ0
export δ0_interval1_deg, δ0_interval2_deg, δ0_deg
export η0_interval1, η0_interval2, η0
include("extension19.jl")

end
