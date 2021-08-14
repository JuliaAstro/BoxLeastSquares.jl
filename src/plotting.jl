using RecipesBase

const OBJ_STRING = Dict(
    :likelihood => "log-likeilhood",
    :snr => "S/N"
)

@recipe function f(bls::BLSPeriodogram)
    seriestype --> :line
    yguide --> OBJ_STRING[bls.objective]
    xguide --> "period"

    periods(bls), power(bls)
end