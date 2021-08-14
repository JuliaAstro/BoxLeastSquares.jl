using BoxLeastSquares
using Distributions
using StableRNGs
using Test

rng = StableRNG(129325)

@testset "BoxLeastSquares.jl" begin
    include("bls.jl")
end
