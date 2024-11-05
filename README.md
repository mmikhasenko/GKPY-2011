# PiPiMadrid

[![Build Status](https://github.com/mmikhasenko/PiPiMadrid.jl/workflows/Test/badge.svg)](https://github.com/mmikhasenko/PiPiMadrid.jl/actions)
[![Test workflow status](https://github.com/mmikhasenko/PiPiMadrid.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/PiPiMadrid.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/mmikhasenko/PiPiMadrid.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/PiPiMadrid.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mmikhasenko/PiPiMadrid.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/PiPiMadrid.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

Julia implementation of the $\pi\pi$ scattering amplitudes from Madrid group.

1. Module `GKPY11` implements S,P, and F waves following the ["Pion-pion scattering amplitude, IV"](https://inspirehep.net/literature/889131) paper.
2. Module `PRR19` codes the S-wave from the ["Global parameterization of ππ scattering up to 2 GeV"](https://inspirehep.net/literature/1747223)

## Installation

```julia
pkg> add https://github.com/mmikhasenko/PiPiMadrid.jl.git
```

## Usage

Function `δ0`, `δ1`, and `δ2` give values of the phase shifts in S, P, and F waves, respectively, as functions of the squared mass of the $\pi\pi$ system.

Signature of each function is
```
δ0(s; pars)
```
where the parameters are set to default values given by `s_wave_pars`, `p_wave_pars`, and `f_wave_pars`.

```julia
using PiPiMadrid
using Plots

theme(:wong2, lab = "", frame = :box, grid = false)

let
    support = (0.3, 2.0)
    plot(xlab = "m(ππ) [GeV]",
        ylab = "phase shift [rad]",
        title = "scattering phase shifts", leg = :topleft)
    plot!(e -> s_wave_phase_shift(GKPY11_default(), e^2), support..., label = "δ0 GKPY11")
    plot!(e -> p_wave_phase_shift(GKPY11_default(), e^2), support..., label = "δ1 GKPY11")
    plot!(e -> f_wave_phase_shift(GKPY11_default(), e^2), support..., label = "δ3 GKPY11")
    plot!(e -> s_wave_phase_shift(PRR19_default(), e^2), 0.3, 2.0, label = "δ0 PRR19")
    vspan!([GKPY11_default().S.m_boundary3, support[2]], alpha = 0.2)
end

savefig("~/Documents/phase_shift.png")
```

![phase_shift](https://github.com/user-attachments/assets/aa358cf9-c117-43dd-8a0a-c1db18400289)

The interval below 1.4GeV is constrained by the GKPY model, the higher range is effective parametrization that fits the data.

## How to Cite

Please cite the original paper, [link to inspire](https://inspirehep.net/literature/889131).
You can refer to Julia implementation with a Zenodo DOI, [link to Zenodo](https://doi.org/FIXME).

