using RecipesBase

@recipe function f(bls::BLSPeriodogram)
    seriestype --> :line
    yguide --> string(bls.objective)
    xguide --> "period"

    periods(bls), power(bls)
end