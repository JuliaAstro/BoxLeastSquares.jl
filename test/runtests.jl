using BoxLeastSquares
using Distributions
using StableRNGs
using Test
using Unitful

rng = StableRNG(129325)

function make_data(N=500)
    t = rand(rng, Uniform(0, 10), N)
    y = ones(N)
    dy = rand(rng, Uniform(0.005, 0.010), N)
    period = 2.0
    t0 = 0.5
    duration = 0.16
    depth = 0.2
    mask = @. abs((t - t0 + 0.5 * period) % period - 0.5 * period) < 0.5 * duration
    y[mask] .-= depth
    y .+= dy .* randn(rng, N)
    return t, y, dy, (;period, t0, duration, depth)
end

@testset "BoxLeastSquares.jl" begin
    include("bls.jl")
    include("plotting.jl")
    include("printing.jl")
end
