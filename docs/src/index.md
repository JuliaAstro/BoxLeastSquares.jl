```@meta
CurrentModule = BoxLeastSquares
```

# BoxLeastSquares.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/JuliaAstro/BoxLeastSquares.jl)
[![Build Status](https://github.com/JuliaAstro/BoxLeastSquares.jl/workflows/CI/badge.svg?branch=main)](https://github.com/JuliaAstro/BoxLeastSquares.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/B/BoxLeastSquares.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/JuliaAstro/BoxLeastSquares.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/BoxLeastSquares.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


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

First, import the package

```jldoctest usage
julia> using BoxLeastSquares
```

you can optionally alias the package name, too

```julia
julia> import BoxLeastSquares as BLS
```

now, load some data. If you don't have an estimate of the y error it will default to 1.

```@meta
DocTestSetup = quote
    using StableRNGs
    rng = StableRNG(3351)
    t = 10 .* rand(rng, 1000);
    yerr = 5e-3 .* (rand(rng, 1000) .+ 1)
    P = 2; t0 = 0.5; dur = 0.16; depth = 0.2;
    mask = @. abs((t - t0 + 0.5P) % P - 0.5P) < 0.5dur;
    y = @. ifelse(mask, 1 - depth, 1);
    y .+= yerr .* randn(rng, 1000)
    load_data() = t, y, yerr
end
```

```jldoctest usage
julia> t, y, yerr = load_data(); # load data somehow
```

The primary interface is through the `BLS` method

```jldoctest usage
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
index: 1633
period: 1.99930396919953
duration: 0.16
t0: 0.5001330656464655
depth: 0.19594118110109113 ± 0.0008688097746093883
snr: 225.52828804117118
log-likelihood: 27396.365214805144
```

to extract the parameters in a convenient named tuple use [`BoxLeastSquares.params`](@ref)

```jldoctest usage
julia> BoxLeastSquares.params(result)
(index = 1633, power = 27396.365214805144, period = 1.99930396919953, duration = 0.16, t0 = 0.5001330656464655, depth = 0.19594118110109113, depth_err = 0.0008688097746093883, snr = 225.52828804117118, loglike = 27396.365214805144)
```

The period grid was automatically determined using [`autoperiod`](@ref), but you can supply your own, too:

```jldoctest usage
julia> periods = exp.(range(log(2) - 0.1, log(2) + 0.1, length=1000));

julia> result_fine = BLS(t, y, yerr; duration=0.12:0.01:0.20, periods=periods)
BLSPeriodogram
==============
input dim: 1000
output dim: 1000
period range: 1.809674836071919 - 2.210341836151295
duration range: 0.12 - 0.2
objective: likelihood

parameters
----------
index: 503
period: 2.001001251543549
duration: 0.168
t0: 0.4961330656464656
depth: 0.19466955969052016 ± 0.0008627202098527317
snr: 225.64622628204188
log-likelihood: 27457.6383039924
```

### Unitful.jl

BoxLeastSquares.jl is fully compatible with `Unitful.jl` (although it is not a dependency of the library). For example

```jldoctest usage
julia> using Unitful

julia> tu = t * u"d";

julia> results_units = BLS(tu, y, yerr; duration=(2:0.1:4)u"hr")
BLSPeriodogram
==============
input dim: 1000
output dim: 3343
period range: 0.3333333333333333 d - 4.988348864592586 d
duration range: 2.0 hr - 4.0 hr
objective: likelihood

parameters
----------
index: 2986
period: 2.0019235780121827 d
duration: 3.8000000000000003 hr
t0: 0.4916330656464656 d
depth: 0.19445716575012517 ± 0.0008692454825826517
snr: 223.70799693127577
log-likelihood: 26953.643422397385
```

### Plotting

[`BoxLeastSquares.BLSPeriodogram`](@ref) has plotting shorthands built right in- by default it will plot the period grid and the computed power

```@example usage
using BoxLeastSquares, Unitful, StableRNGs # hide
rng = StableRNG(3351) # hide
t = 10 .* rand(rng, 1000)u"d" # hide
yerr = 5e-3 .* (rand(rng, 1000) .+ 1) # hide
P = 2u"d"; t0 = 0.5u"d"; dur = 0.16u"d"; depth = 0.2; # hide
mask = @. abs((t - t0 + 0.5P) % P - 0.5P) < 0.5dur # hide
y = @. ifelse(mask, 1.0 - depth, 1.0) # hide
y .+= yerr .* randn(rng, 1000) # hide
results_units = BLS(t, y, yerr; duration=(2:0.1:4)u"hr") # hide
using Plots, UnitfulRecipes

plot(results_units, label="")
```

now let's look at how the transit model compares to the data

```@example usage
pars = BoxLeastSquares.params(results_units)
wrap = 0.5 * pars.period
phases = @. (mod(t - pars.t0 + wrap, pars.period) - wrap) / pars.period
inds = sortperm(phases)
model = BoxLeastSquares.model(results_units)

scatter(phases[inds], y[inds], yerr=yerr[inds],
    label="data", xlabel="phase", xlim=(-0.2, 0.2), leg=:bottomright)
plot!(phases[inds], model[inds], lw=3, label="BLS model")
```

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaAstro/BoxLeastSquares.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaAstro/BoxLeastSquares.jl/discussions) and join or open a new topic. If you're having problems with something, open an [issue](https://github.com/JuliaAstro/BoxLeastSquares.jl/issues).
