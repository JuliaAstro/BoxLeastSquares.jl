module BoxLeastSquares

using Statistics

export BLS, autoperiod

# core BLS algorithm and struct
include("bls.jl")
# transit model
include("model.jl")
# pretty printing
include("printing.jl")

end
