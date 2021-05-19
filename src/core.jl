

struct BLSPeriodogram{TT,VT<:AbstractVector{TT},YT<:AbstractVector,ET<:AbstractVector}
    t::VT
    y::YT
    err::ET
    period
    duration
    power
    method::Symbol
end

period(bls::BLSPeriodogram) = bls.period
duration(bls::BLSPeriodogram) = bls.duration
power(bls::BLSPeriodogram) = bls.power

function Base.show(io::IO, ::MIME"text/plain", bls::BLSPeriodogram)
    best_ind = argmax(power(bls))
    best_pow = power(bls)[best_ind]
    best_per = period(bls)[best_ind]
    method = string(bls.method)
    print(io, "BLSPeriodogram\n", "best period: ", best_per, "\n$method at best period: ", best_pow)
end


function autoperiod(t,
                    duration;
                    minimum_n_transit=3,
                    minimum_period=default_minimum_period(duration),
                    maximum_period=default_maximum_period(t, minimum_n_transit),
                    frequency_factor=1.0)

    lo, hi = extrema(t)
    baseline = hi - lo
    df = frequency_factor * minimum(duration) / baseline^2

    if maximum_period < minimum_period
        minimum_period, maximum_period = maximum_period, minimum_period
    end
    minimum_period ≥ 0 || error("minimum period must be positive")

    minimum_frequency = inv(maximum_period)
    maximum_frequency = inv(minimum_period)

    nf = (maximum_frequency - minimum_frequency) ÷ df
    return @. inv(maximum_frequency - df * (0:nf))
end

default_minimum_period(duration) = 2 * maximum(duration)
function default_maximum_period(t, minimum_n_transit)
    minimum_n_transit > 1 || error("minimum number of transits must be greater than 1")
    lo, hi = extrema(t)
    return (hi - lo) / (minimum_n_transit - 1)
end

function bls(t, y, yerr;
               duration,
               objective=:likelihood,
               oversample=10,
               minimum_n_transit=3,
               minimum_period=default_minimum_period(duration),
               maximum_period=default_maximum_period(t, minimum_n_transit),
               frequency_factor=1.0,
               period=autoperiod(t,
                                 duration;
                                 minimum_n_transit=minimum_n_transit,
                                 minimum_period=minimum_period,
                                 maximum_period=maximum_period,
                                 frequency_factor=frequency_factor))
    invvar = @. inv(yerr^2)
    y_res = y .- median(y)

    power = map(period) do P
        best = -Inf
        half_P = 0.5 * P
        min_t = minimum(t)
        for τ in duration
            # compute phase grid
            dp = τ / oversample
            phase = 0:dp:P

            for t0 in phase
                # find in-transit data points
                mask = @. abs((t - min_t - t0 + half_P) % P - half_P) < 0.5 * τ
                imask = .!mask

                # estimate in-transit flux
                invvar_in = sum(invvar[mask])
                y_in = mapreduce(*, +, y[mask], invvar[mask]) / invvar_in
                # estimate out-of-transit flux
                invvar_out = sum(invvar[imask])
                y_out = mapreduce(*, +, y[imask], invvar[imask]) / invvar_out

                # compute best-fit depth
                if objective === :snr
                    depth = y_out - y_in
                    depth_err = sqrt(inv(invvar_in) + inv(invvar_out))
                    snr = depth / depth_err
                    best = max(best, snr)
                elseif objective === :likelihood
                    loglike = -0.5 * mapreduce((y,iv) -> ((y - y_in)^2 - (y - y_out)^2) * iv, +, y[mask], invvar[mask])
                    best = max(best, loglike)
                else
                    error("invalid objective $objective. Should be one of `:snr` or `:likelihood`")
                end
            end
        end
        return best
    end
    return BLSPeriodogram(t, y, yerr, period, duration, power, objective)
end
