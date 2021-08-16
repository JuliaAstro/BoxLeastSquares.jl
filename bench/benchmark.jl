using BenchmarkTools
using BoxLeastSquares
using CSV
using ProgressLogging
using PyCall
using Random
using Statistics

ENV["OMP_NUM_THREADS"] = 1
@info "OMP NUM THREADS=$(ENV["OMP_NUM_THREADS"])"

apt = pyimport("astropy.timeseries")


# generate data
rng = Random.MersenneTwister(444219)
params = (
    depth=0.1,
    period=3,
    duration=0.2
)
function gendata(N)
    t = 20 .* rand(rng, N)
    y = @. 1 - params.depth * (mod(t, params.period) < params.duration)
    yerr = 0.005 .* (rand(rng, N) .+ 1)
    y .+= yerr .* randn(rng, N)
    return t, y, yerr, params
end

function periodogram_astropy(t, y, yerr, periods)
    model = apt.BoxLeastSquares(t, y, dy=yerr)
    periodogram = model.power(periods, params.duration)
    return model, periodogram
end

function periodogram_boxleastsquares(t, y, yerr, periods)
    return BLS(t, y, yerr; duration=params.duration, periods)
end

data_sizes = round.(Int, 10 .^ range(2, 4, length=11))
period_sizes = round.(Int, 10 .^ range(2, 4, length=3))

results = []
@progress for N_data in data_sizes
    t, y, yerr, params = gendata(N_data)
    @progress for N_per in period_sizes
        log_periods = range(log(params.period) - 1, log(params.period) + 1, length=N_per)
        periods = exp.(log_periods)

        bench_jl = @benchmark periodogram_boxleastsquares($t, $y, $yerr, $periods)
        bench_py = @benchmark periodogram_astropy($t, $y, $yerr, $periods)
        
        res = (N_data, N_per, t_astropy=median(bench_py).time/1e9, t_bls=median(bench_jl).time/1e9)
        @info "" res...
        push!(results, res)
    end
end

CSV.write(joinpath(@__DIR__, "benchmark_results.csv"), results)