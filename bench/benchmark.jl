using BenchmarkTools
using BoxLeastSquares
using PyCall
using ProgressLogging
using Random
using Statistics

ENV["OMP_NUM_THREADS"] = 1
@info "OMP NUM THREADS=$(ENV["OMP_NUM_THREADS"])"

apt = pyimport("astropy.timeseries")

# generate data
rng = Random.MersenneTwister(444219)
function gendata(N)
    t = 20 .* rand(rng, N)
    y = @. 1 - 0.1 * ((t % 3) < 0.2)
    y .+= 0.01 .* randn(rng, N)
    yerr = fill(0.01, N)
    return t, y, yerr
end

t, y, yerr = gendata(1000)

periods = autoperiod(t, 0.2)

function periodogram_astropy(t, y, yerr, periods)
    model = apt.BoxLeastSquares(t, y, dy=yerr)
    periodogram = model.power(periods, 0.2)
    return model, periodogram
end

function periodogram_boxleastsquares(t, y, yerr, periods)
    return BLS(t, y, yerr; duration=0.2, periods)
end

data_sizes = round.(Int, 10 .^ range(2, 4, length=11))
period_sizes = round.(Int, 10 .^ range(2, 4, length=3))

results = []
@progress for N_data in data_sizes
    t, y, yerr = gendata(N_data)
    true_period = 3
    @progress for N_per in period_sizes
        log_periods = range(log(true_period) - 1, log(true_period) + 1, length=N_per)
        periods = exp.(log_periods)

        bench_jl = @benchmark periodogram_boxleastsquares($t, $y, $yerr, $periods)
        bench_py = @benchmark periodogram_astropy($t, $y, $yerr, $periods)
        
        res = (N_data, N_per, t_astropy=median(bench_py), t_bls=median(bench_jl))
        @info "" res...
        push!(results, res)
    end
end

CSV.write(joinpath(@__DIR__, "benchmark_results.csv"), results)