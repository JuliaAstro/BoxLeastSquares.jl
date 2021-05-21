

struct BLSPeriodogram{TT,VT<:AbstractVector{TT},YT<:AbstractVector,ET<:AbstractVector,PT,PWT,DT,TOT}
    t::VT
    y::YT
    err::ET
    period::PT
    power::PWT
    duration::DT
    t0::TOT
    method::Symbol
end

period(bls::BLSPeriodogram) = bls.period
duration(bls::BLSPeriodogram) = bls.duration
power(bls::BLSPeriodogram) = bls.power
t0(bls::BLSPeriodogram) = bls.t0

function Base.show(io::IO, ::MIME"text/plain", bls::BLSPeriodogram)
    best_ind = argmax(power(bls))
    best_pow = power(bls)[best_ind]
    best_per = period(bls)[best_ind]
    best_dur = duration(bls)[best_ind]
    best_t0 = t0(bls)[best_ind]
    method = string(bls.method)
    n = length(bls.t)
    m = length(bls.power)
    print(io, "BLSPeriodogram\n==============\ninput dim: ", n, "\noutput dim: ", m, "\n\nparameters\n----------\n", "period: ", best_per, "\nduration: ", best_dur, "\nt0: ", best_t0, "\n$method: ", best_pow)
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

    nf = 1 + (maximum_frequency - minimum_frequency) ÷ df
    return @. inv(maximum_frequency - df * (0:nf))
end

default_minimum_period(duration) = 2 * maximum(duration)
function default_maximum_period(t, minimum_n_transit)
    minimum_n_transit > 1 || error("minimum number of transits must be greater than 1")
    lo, hi = extrema(t)
    return (hi - lo) / (minimum_n_transit - 1)
end

function BLS(t, y, yerr;
               duration,
               objective=:likelihood,
               oversample=10,
               period=autoperiod(t, duration;))
    invvar = @. inv(yerr^2)
    y_res = y .- median(y)
    min_t = minimum(t)

    # set up arrays
    powers = similar(period, float(eltype(period)))
    durations = similar(period, typeof(duration))
    t0s = similar(period)

    for (idx, P) in pairs(period)
        best = -Inf
        half_P = 0.5 * P
        for τ in duration
            # compute phase grid
            dp = τ / oversample
            phase = 0:dp:P

            for t0 in phase
                # find in-transit data points
                invvar_in = invvar_out = zero(eltype(invvar))
                y_in = y_out = zero(eltype(y_res))
                for i in eachindex(t)
                    transiting = abs((t[i] - min_t - t0 + half_P) % P - half_P) < 0.5 * τ
                    if transiting
                        invvar_in += invvar[i]
                        y_in += y_res[i] * invvar[i]
                    else
                        invvar_out += invvar[i]
                        y_out += y_res[i] * invvar[i]
                    end
                end
                y_in /= invvar_in
                y_out /= invvar_out

                # # estimate in-transit flux
                # invvar_in = sum(invvar[mask])
                # y_in = mapreduce(*, +, y_res[mask], invvar[mask]) / invvar_in
                # # estimate out-of-transit flux
                # invvar_out = sum(invvar[imask])
                # y_out = mapreduce(*, +, y_res[imask], invvar[imask]) / invvar_out

                # compute best-fit depth
                if objective === :snr
                    depth = y_out - y_in
                    depth_err = sqrt(inv(invvar_in) + inv(invvar_out))
                    snr = depth / depth_err
                    if snr > best
                        best = snr
                        powers[idx] = snr
                        durations[idx] = τ
                        t0s[idx] = t0
                    end
                elseif objective === :likelihood
                    loglike = 0.0
                    for i in eachindex(t)
                        transiting = abs((t[i] - min_t - t0 + half_P) % P - half_P) < 0.5 * τ
                        transiting || continue
                        loglike += (y_res[i] - y_in)^2 * invvar[i]
                        loglike -= (y_res[i] - y_out)^2 * invvar[i]
                    end
                    loglike *= -0.5
                    if loglike > best
                        best = loglike
                        powers[idx] = loglike
                        durations[idx] = τ
                        t0s[idx] = t0
                    end
                else
                    error("invalid objective $objective. Should be one of `:snr` or `:likelihood`")
                end
            end
        end
    end
    return BLSPeriodogram(t, y, yerr, period, powers, durations, t0s, objective)
end
