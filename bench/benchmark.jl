using BenchmarkTools
using BoxLeastSquares
using BoxLeastSquares: period, power
using PyCall
using Random

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

function periodogram_astropy(t, y, yerr)
    model = apt.BoxLeastSquares(t, y, dy=yerr)
    periodogram = model.autopower(0.2)
    return periodogram["period"], periodogram["power"]
end

function periodogram_boxleastsquares(t, y, yerr)
    periodogram = BLS(t, y, yerr; duration=0.2)
    return period(periodogram), power(periodogram)
end

@benchmark periodogram_astropy($t, $y, $yerr)
@benchmark periodogram_boxleastsquares($t, $y, $yerr)
