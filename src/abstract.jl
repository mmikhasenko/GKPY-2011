abstract type PiPiBase end

"""
    s_wave_amplitude(model::PiPiBase, s)

Returns the S-wave scattering amplitude for a given model at energy squared `s` (in GeV²).
Subtypes of `PiPiBase` should implement this function to provide a model-specific calculation
of the S-wave amplitude.
"""
function s_wave_amplitude(model::PiPiBase, s; pars = model.S)
    throw(MethodError(s_wave_amplitude, (model, s)))
end

"""
    s_wave_phase_shift(model::PiPiBase, s)

Returns the S-wave phase shift in radians for a given model at energy squared `s` (in GeV²).
This function should be implemented by subtypes to provide model-specific calculations.
"""
function s_wave_phase_shift(model::PiPiBase, s; pars = model.S)
    throw(MethodError(s_wave_phase_shift, (model, s)))
end

"""
    s_wave_elasticity(model::PiPiBase, s)

Returns the S-wave elasticity for a given model at energy squared `s` (in GeV²).
Elasticity values indicate the fraction of the total cross-section due to elastic scattering.
Subtypes should implement this function for model-specific elasticity calculations.
"""
function s_wave_elasticity(model::PiPiBase, s; pars = model.S)
    throw(MethodError(s_wave_elasticity, (model, s)))
end

"""
    s_wave_coupled_channel_matrix(model::PiPiBase, s)

Returns the coupled channel matrix for the S-wave for a given model at energy squared `s` (in GeV²).
This function is used for models with coupled-channel interactions.
"""
function s_wave_coupled_channel_matrix(model::PiPiBase, s; pars = model.S)
    throw(MethodError(s_wave_coupled_channel_matrix, (model, s)))
end

"""
    p_wave_amplitude(model::PiPiBase, s)

Returns the P-wave scattering amplitude for a given model at energy squared `s` (in GeV²).
Subtypes of `PiPiBase` should implement this function to provide a model-specific calculation
of the P-wave amplitude.
"""
function p_wave_amplitude(model::PiPiBase, s; pars = model.P)
    throw(MethodError(p_wave_amplitude, (model, s)))
end

"""
    p_wave_phase_shift(model::PiPiBase, s)

Returns the P-wave phase shift in radians for a given model at energy squared `s` (in GeV²).
This function should be implemented by subtypes to provide model-specific calculations.
"""
function p_wave_phase_shift(model::PiPiBase, s; pars = model.P)
    throw(MethodError(p_wave_phase_shift, (model, s)))
end

"""
    p_wave_elasticity(model::PiPiBase, s)

Returns the P-wave elasticity for a given model at energy squared `s` (in GeV²).
Elasticity values indicate the fraction of the total cross-section due to elastic scattering.
Subtypes should implement this function for model-specific elasticity calculations.
"""
function p_wave_elasticity(model::PiPiBase, s; pars = model.P)
    throw(MethodError(p_wave_elasticity, (model, s)))
end

"""
    d_wave_amplitude(model::PiPiBase, s)

Returns the D-wave scattering amplitude for a given model at energy squared `s` (in GeV²).
This function should be implemented by subtypes to provide a model-specific D-wave amplitude.
"""
function d_wave_amplitude(model::PiPiBase, s; pars = model.D)
    throw(MethodError(d_wave_amplitude, (model, s)))
end

"""
    d_wave_phase_shift(model::PiPiBase, s)

Returns the F-wave phase shift in radians for a given model at energy squared `s` (in GeV²).
This function should be implemented by subtypes to provide model-specific calculations.
"""
function d_wave_phase_shift(model::PiPiBase, s; pars = model.F)
    throw(MethodError(d_wave_phase_shift, (model, s)))
end

"""
    f_wave_amplitude(model::PiPiBase, s)

Returns the F-wave scattering amplitude for a given model at energy squared `s` (in GeV²).
This function should be implemented by subtypes to provide a model-specific F-wave amplitude.
"""
function f_wave_amplitude(model::PiPiBase, s; pars = model.F)
    throw(MethodError(f_wave_amplitude, (model, s)))
end

"""
    f_wave_phase_shift(model::PiPiBase, s)

Returns the F-wave phase shift in radians for a given model at energy squared `s` (in GeV²).
This function should be implemented by subtypes to provide model-specific calculations.
"""
function f_wave_phase_shift(model::PiPiBase, s; pars = model.F)
    throw(MethodError(f_wave_phase_shift, (model, s)))
end
