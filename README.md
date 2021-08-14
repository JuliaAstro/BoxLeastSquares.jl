# BoxLeastSquares.jl

[![Build Status](https://github.com/JuliaAstro/BoxLeastSquares.jl/workflows/CI/badge.svg?branch=main)](https://github.com/JuliaAstro/BoxLeastSquares.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/B/BoxLeastSquares.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/JuliaAstro/BoxLeastSquares.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/BoxLeastSquares.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaAstro.github.io/BoxLeastSquares.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaAstro.github.io/BoxLeastSquares.jl/dev)

Box-least-squares (BLS) periodograms in pure Julia. 


## Installation

To install use [Pkg](https://julialang.github.io/Pkg.jl/v1/managing-packages/). From the REPL, press `]` to enter Pkg-mode

```julia
pkg> add BoxLeastSquares
```
If you want to use the most up-to-date version of the code, check it out from `main`

```julia
pkg> add BoxLeastSquares#main
```

## Usage

First, import the package, optionally aliasing the package name

```julia
julia> using BoxLeastSquares

julia> import BoxLeastSquares as BLS
```

now, load some data. If you don't have an estimate of the y error it will default to 1.

```julia
julia> t, y, yerr = # produce data
```

The primary interface is through the `BLS` method

```julia
julia> result = BLS(t, y, yerr; duration=0.16)
BLSPeriodogram
==============
input dim: 1000
output dim: 1820
period range: 0.32 - 5.014724142709022
duration range: 0.16 - 0.16
objective: likelihood

parameters
----------
period: 1.99930396919953
duration: 0.16
t0: 0.5001330656464655
depth: 0.19594118110109113 ± 0.0008688097746093883
snr: 225.52828804117118
log-likelihood: 27396.365214805144
```

The transit parameters at the maximum power can be retrieved using `BoxLeastSquares.params`

```julia
julia> BoxLeastSquares.params(result)
(power = 27396.365214805144, period = 1.99930396919953, duration = 0.16, t0 = 0.5001330656464655, depth = 0.19594118110109113, depth_err = 0.0008688097746093883, snr = 225.52828804117118, loglike = 27396.365214805144)
```

The period grid was automatically determined using `autoperiod`, but you can supply your own, too:

```julia
julia> periods = range(log(2) - 0.1, log(2) + 0.1, length=1000);

julia> result_fine = BLS(t, y, yerr; duration=0.12:0.20, periods=periods)
BLSPeriodogram
==============
input dim: 1000
output dim: 1000
period range: 0.5931471805599453 - 0.7931471805599453
duration range: 0.12 - 0.12
objective: likelihood

parameters
----------
period: 0.668822856235621
duration: 0.12
t0: 0.49613306564646553
depth: 0.06949975333493467 ± 0.0005991935962031148
snr: 115.98881192210811
log-likelihood: 8087.810929044984
```

### Unitful.jl

BoxLeastSquares.jl is fully compatible with the `Unitful.jl` (although it is not a dependency of the library). For example

```julia
julia> using Unitful

julia> tu = t * u"d";

julia> results_units = BLS(tu, y, yerr; duration=(10:14)u"hr")
BLSPeriodogram
==============
input dim: 1000
output dim: 158
period range: 1.1666666666666665 d - 5.010959571954422 d
duration range: 10 hr - 14 hr
objective: likelihood

parameters
----------
period: 2.04681966153818 d
duration: 10 hr
t0: 0.3791330656464655 d
depth: 0.0788086190813451 ± 0.0005822150269728557
snr: 135.35998802899218
log-likelihood: 11186.069545672593
```

## Contributing and Support
 
If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaAstro/BoxLeastSquares.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaAstro/BoxLeastSquares.jl/discussions) and join or open a new topic. If you're having problems with something, open an [issue](https://github.com/JuliaAstro/BoxLeastSquares.jl/issues).
