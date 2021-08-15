
@testset "BLSPeriodogram printing" begin
    t, y, dy, params = make_data()
    result = BLS(t, y, dy; params.duration)

    N_in = length(t)
    N_out = length(result.periods)
    T = eltype(t)
    V = eltype(y)
    # shortform
    @test sprint(show, result) == "$N_out-element BLSPeriodogram{$T,$V}"
    
    # longform
    minP, maxP = extrema(result.periods)
    pars = BoxLeastSquares.params(result)
    @test repr("text/plain", result) == """
    BLSPeriodogram
    ==============
    input dim: $N_in
    output dim: $N_out
    period range: $minP - $maxP
    duration range: $(params.duration) - $(params.duration)
    objective: $(string(result.objective))

    parameters
    ----------
    index: $(pars.index)
    period: $(pars.period)
    duration: $(pars.duration)
    t0: $(pars.t0)
    depth: $(pars.depth) ± $(pars.depth_err)
    snr: $(pars.snr)
    log-likelihood: $(pars.loglike)"""
end