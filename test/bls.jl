using Unitful: unit


@testset "autoperiod self-consistency" begin
    t, y, dy, params = make_data()
    durations = params.duration .+ range(-0.1, 0.1, length=3)

    period = @inferred autoperiod(t, durations)
    results1 = @inferred BLS(t, y, dy; duration=durations, periods=period)
    results2 = @inferred BLS(t, y, dy; duration=durations)
    @test BoxLeastSquares.power(results1) ≈ BoxLeastSquares.power(results2)

end

@testset "model correctness ($obj)" for obj in [:likelihood, :snr]
    t, y, dy, params = make_data(1000)

    periods = exp.(range(log(params.period) - 0.1, log(params.period) + 0.1, length=1000))

    results = @inferred BLS(t, y, dy; params.duration, periods, objective=obj)
    best_params = BoxLeastSquares.params(results)

    @test best_params.period ≈ params.period atol = 0.01
    @test best_params.t0 ≈ params.t0 atol = 0.01
    @test best_params.duration ≈ params.duration atol = 0.01
    @test best_params.depth ≈ params.depth atol = best_params.depth_err
end

@testset "transit model" begin
    t, y, dy, params = make_data()

    # compute model using linear regression
    A = zeros(length(t), 2)
    dt = @. abs((t - params.t0 + 0.5 * params.period) % params.period - 0.5 * params.period)
    intransit = @. dt < 0.5 * params.duration
    A[.!intransit, 1] .= 1
    A[intransit, 2] .= 1
    w = (A' * (A ./ dy .^2)) \ (A' * (y ./ dy .^2))
    model_true = A * w

    # model = @inferred BoxLeastSquares.model(t, y, dy; params...)
    model = BoxLeastSquares.model(t, y, dy; params...)
    @test model ≈ model_true
end

@testset "Unitful.jl integration" begin
    @testset "time units" begin
        t, y, dy, params = make_data()

        tu = t * u"d"
        duration = 0.16u"d" |> u"hr"

        results = BLS(tu, y, dy; duration)
        @test unit(eltype(results.periods)) == u"d"
        @test unit(eltype(results.t0)) == u"d"
        @test unit(eltype(results.duration)) == u"hr"
        @test unit(eltype(results.duration_in)) == u"hr"
    end
    @testset "flux units" begin
        t, y, dy, params = make_data()

        yu = y * u"J"
        dyu = dy * u"J"

        results = BLS(t, yu, dyu; params.duration)
        @test unit(eltype(results.depth)) == u"J"
    end
    @testset "flux and time units" begin
        t, y, dy, params = make_data()

        tu = t * u"d"
        yu = y * u"J"
        dyu = dy * u"J"
        duration = 0.16u"d" |> u"hr"

        results = BLS(tu, yu, dyu; duration)
        @test unit(eltype(results.periods)) == u"d"
        @test unit(eltype(results.t0)) == u"d"
        @test unit(eltype(results.duration)) == u"hr"
        @test unit(eltype(results.duration_in)) == u"hr"
        @test unit(eltype(results.depth)) == u"J"
    end
end
