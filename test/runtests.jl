using AstroTime
using BoxLeastSquares
using Dates
using Distributions
using StableRNGs
using Test
using Unitful

rng = StableRNG(129325)

@testset "BoxLeastSquares.jl" begin  include("bls.jl") end
