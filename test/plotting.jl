using RecipesBase: apply_recipe

@testset "BLSPeriodogram plotting" begin
    t, y, dy, params = make_data()
    result = BLS(t, y, dy; params.duration)
    recipes = apply_recipe(Dict{Symbol,Any}(), result)
    for rec in recipes
        @test getfield(rec, 1) == Dict{Symbol,Any}(
            :seriestype => :line,
            :xguide => "period",
            :yguide => "log-likeilhood")

            # test to make sure our position is correct (should be +0.5 the given (x,y))
        period = rec.args[1]
        power = rec.args[2]
        
        @test period == result.periods
        @test power == result.power
    end
end