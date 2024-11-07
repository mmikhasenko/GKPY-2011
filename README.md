# PiPiMadrid

[![Build Status](https://github.com/mmikhasenko/PiPiMadrid.jl/workflows/Test/badge.svg)](https://github.com/mmikhasenko/PiPiMadrid.jl/actions)
[![Test workflow status](https://github.com/mmikhasenko/PiPiMadrid.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/PiPiMadrid.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/mmikhasenko/PiPiMadrid.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/PiPiMadrid.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mmikhasenko/PiPiMadrid.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/PiPiMadrid.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

`PiPiMadrid.jl` is a Julia package providing implementations of two-particle scattering amplitudes for pion-pion interactions based on models developed by the Madrid Group.
The package includes methods to access these amplitudes in different partial waves.

Two models are available:

1. Module `GKPY11` implements S,P, and F waves following the ["Pion-pion scattering amplitude, IV"](https://inspirehep.net/literature/889131) paper.
2. Module `PRR19` codes the S-wave from the ["Global parameterization of ππ scattering up to 2 GeV"](https://inspirehep.net/literature/1747223)

## Functionality

The package exports the following main functions for each wave channel:

- **Amplitude**: `s_wave_amplitude`, `p_wave_amplitude`, `d_wave_amplitude`, `f_wave_amplitude`
- **Phase Shift**: `s_wave_phase_shift`, `p_wave_phase_shift`, `d_wave_phase_shift`, `f_wave_phase_shift`
- **Elasticity**: `s_wave_elasticity`, `p_wave_elasticity`

All functions have the same signature: `property(model, s; pars)`, where a type of the `model` parameter is used for dispatch, the default value of the `pars` is set to a component of the model, like `model.S` or `model.P`.
The `s` is a Mandelstam variable, `s = m^2`, where `m` is a mass of the pion pair.

Not every method is defined for each model. See a section on contributing to add more methods.

## Installation

Install the package using Julia’s package manager:

```julia
] add https://github.com/mmikhasenko/PiPiMadrid.jl.git
```

## Usage

After installation, `PiPiMadrid.jl` provides functions to evaluate pion-pion scattering amplitudes.
The following example demonstrates plotting the S-wave amplitude using the `GKPY11` model.

```julia
using PiPiMadrid
using Plots

theme(:wong2, lab = "", frame = :box, grid = false)

let
    support = (0.3, 2.0)
    plot(xlab = "m(ππ) [GeV]",
        ylab = "phase shift [rad]",
        title = "scattering phase shifts", leg = :topleft)
    plot!(e -> s_wave_phase_shift(GKPY11_default, e^2), support..., label = "δ0 GKPY11")
    plot!(e -> p_wave_phase_shift(GKPY11_default, e^2), support..., label = "δ1 GKPY11")
    plot!(e -> f_wave_phase_shift(GKPY11_default, e^2), support..., label = "δ3 GKPY11")
    plot!(e -> s_wave_phase_shift(PRR19_default, e^2), 0.3, 2.0, label = "δ0 PRR19")
    vspan!([GKPY11_default.S.m_boundary3, support[2]], alpha = 0.2)
end
```
![phase_shift](https://github.com/user-attachments/assets/aa358cf9-c117-43dd-8a0a-c1db18400289)

Elasticity is computed using `s_wave_elasticity` function.

The interval below 1.4GeV is constrained by forward Roy-Steiner equations, the higher range is effective parametrization that fits the data.

## Contributing

Contributions are welcome to extend `PiPiMadrid.jl`'s functionality, including additional wave channels or improvements.

1. **Issues**: If you find bugs, missing features, or have suggestions, please open an issue.
2. **Pull Requests**: For feature additions or fixes, fork the repository, implement your changes in a branch, and submit a pull request. Please document any new functions or parameters.

## How to Cite

Please cite the original paper, [link to inspire](https://inspirehep.net/literature/889131).
You can refer to Julia implementation with a Zenodo DOI, [link to Zenodo](https://doi.org/FIXME).

