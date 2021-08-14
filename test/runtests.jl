using BoxLeastSquares
using Distributions
using StableRNGs
using Test

rng = StableRNG(8462852)

@testset "BoxLeastSquares.jl" begin
    include("bls.jl")
end
