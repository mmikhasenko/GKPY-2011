# GKPY11

[![Build Status](https://github.com/mmikhasenko/GKPY11.jl/workflows/Test/badge.svg)](https://github.com/mmikhasenko/GKPY11.jl/actions)
[![Test workflow status](https://github.com/mmikhasenko/GKPY11.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/GKPY11.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/mmikhasenko/GKPY11.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/GKPY11.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mmikhasenko/GKPY11.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/GKPY11.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

Julia implementation of the $\pi\pi$ scattering amplitudes in S,P, and F waves following the GKPY11 model, from the ["Pion-pion scattering amplitude, IV"](https://inspirehep.net/literature/889131).


## Installation

```julia
pkg> add https://github.com/mmikhasenko/GKPY11.jl.git
```

## Usage

Function `δ0`, `δ1`, and `δ2` give values of the phase shifts in S, P, and F waves, respectively, as functions of the squared mass of the $\pi\pi$ system.

Signature of each function is
```
δ0(s; pars)
```
where the parameters are set to defualt values given by `s_wave_pars`, `p_wave_pars`, and `f_wave_pars`.

```julia
using GKPY11
using Plots

let
    support = (0.1, 2.0)
    plot(xlab = "m(ππ) [GeV]",
        ylab = "phase shift [rad]",
        title = "scattering phase shifts", leg=:topleft)
    plot!(e -> δ0(e^2), support..., label = "δ0", lw = 2)
    plot!(e -> δ1(e^2), support..., label = "δ1", lw = 2)
    plot!(e -> δ3(e^2), support..., label = "δ3", lw = 2)
    vspan!([GKPY11.s_wave_pars.m_boundary3, support[2]], alpha=0.1)
end
```


## How to Cite

Please cite the original paper, [link to inspire](https://inspirehep.net/literature/889131).
You can refer to Julia implementation with a Zenodo DOI, [link to Zenodo](https://doi.org/FIXME).

