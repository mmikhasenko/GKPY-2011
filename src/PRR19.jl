module PRR19
using Parameters

import ..iϵ
import ..mπ, ..mK, ..mη
import ..breakup, ..σ, ..J_bar
import ..conformal_w1
import ..GKPY11: cotδ0_interval1

export CheybyshevT
include("cheybyshev.jl")

export t0_interval1, t0_interval2, t0
export δ0_interval1, δ0_interval2, δ0
export δ0_interval1_deg, δ0_interval2_deg, δ0_deg
export η0_interval1, η0_interval2, η0
include("extension19.jl")

end
