module PiPiMadrid

using Parameters
using ForwardDiff: gradient
include("kinematics.jl")

export s_wave_amplitude, s_wave_phase_shift, s_wave_elasticity, s_wave_coupled_channel_matrix
export p_wave_amplitude, p_wave_phase_shift, p_wave_elasticity
export d_wave_amplitude, d_wave_phase_shift
export f_wave_amplitude, f_wave_phase_shift
include("abstract.jl")

# GKPY11 model
export GKPY11_default
include("GKPY11.jl")
include("odd-waves.jl")
include("s-wave.jl")

# PRR19 model
export PRR19_default
include("PRR19.jl")
include("cheybyshev.jl")
include("extension19.jl")


end
